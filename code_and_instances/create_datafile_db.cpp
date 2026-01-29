#include <windows.h>
#include <sqlext.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <cstdio>

// g++ -g "create_datafile_db.cpp" -o "create_datafile_db.exe" -lodbc32 -lodbccp32


struct sql_params{
    SQLHENV env;
    SQLHDBC dbc;
    SQLHSTMT stmt; 
};


// ------------------- TXT file reading functions --------------------
static std::string odbc_diag(SQLSMALLINT handleType, SQLHANDLE handle) {
    SQLCHAR state[6]{};
    SQLCHAR msg[1024]{};
    SQLINTEGER nativeErr = 0;
    SQLSMALLINT msgLen = 0;

    if (SQLGetDiagRecA(handleType, handle, 1, state, &nativeErr, msg, sizeof(msg), &msgLen) == SQL_SUCCESS) {
        return std::string((char*)state) + " | " + std::string((char*)msg);
    }
    return "No diagnostics";
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> out;
    out.reserve(64);
    size_t start = 0;
    while (true) {
        size_t pos = line.find('\t', start);
        if (pos == std::string::npos) { out.push_back(line.substr(start)); break; }
        out.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    // drop trailing empty tokens caused by line ending with tabs
    while (!out.empty() && out.back().empty()) out.pop_back();
    return out;
}

static SQLSMALLINT to_bit(const std::string& s) {
    // accepts "0/1", "true/false" etc if needed
    return (s == "1" || s == "true" || s == "TRUE") ? 1 : 0;
}

static long to_long(const std::string& s) {
    return s.empty() ? 0L : std::stol(s);
}

static double to_double(const std::string& s) {
    return s.empty() ? 0.0 : std::stod(s);
}

static void check_rc(SQLRETURN rc, const char* where) {
    if (rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO) return;
    throw std::runtime_error(std::string("ODBC error at: ") + where);
}

static void check_rc_diag(SQLRETURN rc, const char* where, SQLSMALLINT handleType, SQLHANDLE handle) {
    if (rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO) return;
    throw std::runtime_error(std::string(where) + " | " + odbc_diag(handleType, handle));
}

// ---------------- ODBC Functions -------------------------
std::string produce_conn_str(const std::string& server, const std::string& database) {
    // Note: Encrypt=no for local dev; if you want Encrypt=yes, add TrustServerCertificate=yes
    return "Driver={ODBC Driver 18 for SQL Server};"
           "Server=" + server + ";"
           "Database=" + database + ";"
           "Trusted_Connection=yes;"
           "Encrypt=no;";
};


sql_params connect_to_db(std::string connStr, sql_params sql_param){
    SQLHENV env = sql_param.env;
    SQLHDBC dbc = sql_param.dbc;
    SQLHSTMT stmt = sql_param.stmt;

    // Allocate environment
    SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
    SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);

    // Allocate connection
    SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);

    // Connect using connection string
    SQLRETURN rc = SQLDriverConnectA(
        dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );
    check_rc(rc, "SQLDriverConnectA");

    SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);

    return {env, dbc, stmt};
};


sql_params cleanup_sql_query(sql_params sql_param){
    SQLHENV env = sql_param.env;
    SQLHDBC dbc = sql_param.dbc;
    SQLHSTMT stmt = sql_param.stmt;

    SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);

    return {env, dbc, stmt};
};


sql_params free_up_statement_for_new_query_type(sql_params sql_param){
    SQLHENV env = sql_param.env;
    SQLHDBC dbc = sql_param.dbc;
    SQLHSTMT stmt = sql_param.stmt;

    SQLFreeStmt(stmt, SQL_CLOSE);        // close any open cursor/results
    SQLFreeStmt(stmt, SQL_UNBIND);       // unbind result columns
    SQLFreeStmt(stmt, SQL_RESET_PARAMS); // unbind parameters

    return {env, dbc, stmt};
};


// 
void create_results_table_if_missing(const sql_params& sql_param)
{
    SQLHSTMT stmt = sql_param.stmt;

    const std::string sql = R"SQL(
IF OBJECT_ID(N'dbo.results_summary', N'U') IS NULL
BEGIN
    CREATE TABLE dbo.results_summary (
        id INT IDENTITY(1,1) PRIMARY KEY,
        instance_id INT NOT NULL,
        file_name NVARCHAR(255) NOT NULL,

        transfers BIT,
        relative_mrt BIT,

        mrt_factor FLOAT,
        speed_factor_pt FLOAT,
        pt_interval FLOAT,

        iterations INT,
        iter_without_impr INT,

        max_removal_pct FLOAT,
        max_deterior_pct FLOAT,

        clustered_removal BIT,
        idarp_adapted BIT,
        without_operator BIT,

        terminals_per_node INT,

        priority_non_pt FLOAT,
        priority_pt FLOAT,
        priority_long_dist FLOAT,
        priority_short_dist FLOAT,

        current_iteration INT,
        total_distance FLOAT,
        computation_time FLOAT,

        avg_actual_ride_time FLOAT,
        avg_actual_ride_time_PT FLOAT,
        avg_actual_ride_time_noPT FLOAT,

        avg_excess_ride_time FLOAT,
        avg_excess_ride_time_PT FLOAT,
        avg_excess_ride_time_noPT FLOAT,

        avg_detour_factor FLOAT,
        avg_detour_factor_PT FLOAT,
        avg_detour_factor_noPT FLOAT,

        avg_waiting_time_PT FLOAT,
        pt_usage_ratio FLOAT,

        accept_random_removal INT,
        accept_worst_removal INT,
        accept_related_removal INT,
        accept_route_removal INT,

        accept_random_order_insertion INT,
        accept_greedy_insertion INT,
        accept_two_regret_insertion INT,

        improve_random_removal INT,
        improve_worst_removal INT,
        improve_related_removal INT,
        improve_route_removal INT,

        improve_random_order_insertion INT,
        improve_greedy_insertion INT,
        improve_two_regret_insertion INT
    );
END
)SQL";

SQLRETURN rc = SQLExecDirectA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
check_rc_diag(rc, "SQLExecDirectA(create table)", SQL_HANDLE_STMT, stmt);

};


// Import TSV into dbo.results_summary (id is IDENTITY, so not inserted).
void import_results_summary_tsv(const std::string& path, int instance_id, sql_params& sql_param)
{
    SQLHSTMT stmt = sql_param.stmt;

    // Reset stmt if you reuse it
    SQLFreeStmt(stmt, SQL_CLOSE);
    SQLFreeStmt(stmt, SQL_UNBIND);
    SQLFreeStmt(stmt, SQL_RESET_PARAMS);

    // Start transaction (per file)
    check_rc(SQLExecDirectA(stmt, (SQLCHAR*)"BEGIN TRANSACTION", SQL_NTS), "BEGIN TRANSACTION");

    try {

    // NOTE: columns list must match your dbo.results_summary schema
    const std::string sql = R"SQL(
INSERT INTO dbo.results_summary (
    instance_id, file_name,
    transfers, relative_mrt,
    mrt_factor, speed_factor_pt, pt_interval,
    iterations, iter_without_impr,
    max_removal_pct, max_deterior_pct,
    clustered_removal, idarp_adapted, without_operator,
    terminals_per_node,
    priority_non_pt, priority_pt, priority_long_dist, priority_short_dist,
    current_iteration, total_distance, computation_time,
    avg_actual_ride_time, avg_actual_ride_time_PT, avg_actual_ride_time_noPT,
    avg_excess_ride_time, avg_excess_ride_time_PT, avg_excess_ride_time_noPT,
    avg_detour_factor, avg_detour_factor_PT, avg_detour_factor_noPT,
    avg_waiting_time_PT, pt_usage_ratio,
    accept_random_removal, accept_worst_removal, accept_related_removal, accept_route_removal,
    accept_random_order_insertion, accept_greedy_insertion, accept_two_regret_insertion,
    improve_random_removal, improve_worst_removal, improve_related_removal, improve_route_removal,
    improve_random_order_insertion, improve_greedy_insertion, improve_two_regret_insertion
) VALUES (
    ?, ?,  ?, ?,  ?, ?, ?,  ?, ?,  ?, ?,  ?, ?, ?,  ?,  ?, ?, ?, ?,  ?, ?, ?,  ?, ?, ?,  ?, ?, ?,  ?, ?, ?,  ?, ?, 
    ?, ?, ?, ?,  ?, ?, ?,  ?, ?, ?, ?,  ?, ?, ?
);
)SQL";

    check_rc(SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS), "SQLPrepareA");

    // ---- bind variables (one set reused for every row) ----
    SQLINTEGER p_instance_id = instance_id;

    char file_buf[256]{};
    SQLLEN  file_len = SQL_NTS;

    SQLSMALLINT transfers=0, relative_mrt=0, clustered_removal=0, idarp_adapted=0, without_operator=0;

    double mrt_factor=0, speed_factor_pt=0, pt_interval=0;
    SQLINTEGER iterations=0, iter_without_impr=0;
    double max_removal_pct=0, max_deterior_pct=0;
    SQLINTEGER terminals_per_node=0;

    double priority_non_pt=0, priority_pt=0, priority_long_dist=0, priority_short_dist=0;

    SQLINTEGER current_iteration=0;
    double total_distance=0, computation_time=0;

    double avg_actual_ride_time=0, avg_actual_ride_time_PT=0, avg_actual_ride_time_noPT=0;
    double avg_excess_ride_time=0, avg_excess_ride_time_PT=0, avg_excess_ride_time_noPT=0;
    double avg_detour_factor=0, avg_detour_factor_PT=0, avg_detour_factor_noPT=0;
    double avg_waiting_time_PT=0, pt_usage_ratio=0;

    SQLINTEGER accept_random_removal=0, accept_worst_removal=0, accept_related_removal=0, accept_route_removal=0;
    SQLINTEGER accept_random_order_insertion=0, accept_greedy_insertion=0, accept_two_regret_insertion=0;

    SQLINTEGER improve_random_removal=0, improve_worst_removal=0, improve_related_removal=0, improve_route_removal=0;
    SQLINTEGER improve_random_order_insertion=0, improve_greedy_insertion=0, improve_two_regret_insertion=0;

    int k = 1;
    check_rc(SQLBindParameter(stmt, k++, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0,0, &p_instance_id, 0, nullptr), "Bind instance_id");
    check_rc(SQLBindParameter(stmt, k++, SQL_PARAM_INPUT, SQL_C_CHAR,  SQL_VARCHAR, sizeof(file_buf),0, file_buf, sizeof(file_buf), &file_len), "Bind file");

    // BIT flags as SQL_C_SSHORT -> SQL_BIT works fine
    auto bind_bit = [&](SQLSMALLINT* v, const char* where){
        check_rc(SQLBindParameter(stmt, k++, SQL_PARAM_INPUT, SQL_C_SSHORT, SQL_BIT, 0,0, v, 0, nullptr), where);
    };
    auto bind_int = [&](SQLINTEGER* v, const char* where){
        check_rc(SQLBindParameter(stmt, k++, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0,0, v, 0, nullptr), where);
    };
    auto bind_dbl = [&](double* v, const char* where){
        check_rc(SQLBindParameter(stmt, k++, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_DOUBLE, 0,0, v, 0, nullptr), where);
    };

    bind_bit(&transfers, "Bind transfers");
    bind_bit(&relative_mrt, "Bind relative_mrt");

    bind_dbl(&mrt_factor, "Bind mrt_factor");
    bind_dbl(&speed_factor_pt, "Bind speed_factor_pt");
    bind_dbl(&pt_interval, "Bind pt_interval");

    bind_int(&iterations, "Bind iterations");
    bind_int(&iter_without_impr, "Bind iter_without_impr");

    bind_dbl(&max_removal_pct, "Bind max_removal_pct");
    bind_dbl(&max_deterior_pct, "Bind max_deterior_pct");

    bind_bit(&clustered_removal, "Bind clustered_removal");
    bind_bit(&idarp_adapted, "Bind idarp_adapted");
    bind_bit(&without_operator, "Bind without_operator");

    bind_int(&terminals_per_node, "Bind terminals_per_node");

    bind_dbl(&priority_non_pt, "Bind priority_non_pt");
    bind_dbl(&priority_pt, "Bind priority_pt");
    bind_dbl(&priority_long_dist, "Bind priority_long_dist");
    bind_dbl(&priority_short_dist, "Bind priority_short_dist");

    bind_int(&current_iteration, "Bind current_iteration");
    bind_dbl(&total_distance, "Bind total_distance");
    bind_dbl(&computation_time, "Bind computation_time");

    bind_dbl(&avg_actual_ride_time, "Bind avg_actual_ride_time");
    bind_dbl(&avg_actual_ride_time_PT, "Bind avg_actual_ride_time_PT");
    bind_dbl(&avg_actual_ride_time_noPT, "Bind avg_actual_ride_time_noPT");

    bind_dbl(&avg_excess_ride_time, "Bind avg_excess_ride_time");
    bind_dbl(&avg_excess_ride_time_PT, "Bind avg_excess_ride_time_PT");
    bind_dbl(&avg_excess_ride_time_noPT, "Bind avg_excess_ride_time_noPT");

    bind_dbl(&avg_detour_factor, "Bind avg_detour_factor");
    bind_dbl(&avg_detour_factor_PT, "Bind avg_detour_factor_PT");
    bind_dbl(&avg_detour_factor_noPT, "Bind avg_detour_factor_noPT");

    bind_dbl(&avg_waiting_time_PT, "Bind avg_waiting_time_PT");
    bind_dbl(&pt_usage_ratio, "Bind pt_usage_ratio");

    bind_int(&accept_random_removal, "Bind accept_random_removal");
    bind_int(&accept_worst_removal, "Bind accept_worst_removal");
    bind_int(&accept_related_removal, "Bind accept_related_removal");
    bind_int(&accept_route_removal, "Bind accept_route_removal");

    bind_int(&accept_random_order_insertion, "Bind accept_random_order_insertion");
    bind_int(&accept_greedy_insertion, "Bind accept_greedy_insertion");
    bind_int(&accept_two_regret_insertion, "Bind accept_two_regret_insertion");

    bind_int(&improve_random_removal, "Bind improve_random_removal");
    bind_int(&improve_worst_removal, "Bind improve_worst_removal");
    bind_int(&improve_related_removal, "Bind improve_related_removal");
    bind_int(&improve_route_removal, "Bind improve_route_removal");

    bind_int(&improve_random_order_insertion, "Bind improve_random_order_insertion");
    bind_int(&improve_greedy_insertion, "Bind improve_greedy_insertion");
    bind_int(&improve_two_regret_insertion, "Bind improve_two_regret_insertion");

    // ---- read file & execute per line ----
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty()) continue;

        auto t = split_tabs(line);

        // Expected tokens: 1 (file) + 44 numeric = 45 total
        // If your lines differ, print t.size() and adjust.

        std::cout << "Row size is :" << t.size() << "\n";

        if (t.size() != 126) {
            throw std::runtime_error("Bad row (different number of columns) at line " + std::to_string(line_no) +
                                     " got " + std::to_string(t.size()));
        }

        // 0: file
        std::snprintf(file_buf, sizeof(file_buf), "%s", t[0].c_str());
        file_len = SQL_NTS;

        int i = 1;
        transfers        = to_bit(t[i++]);
        relative_mrt     = to_bit(t[i++]);
        mrt_factor       = to_double(t[i++]);
        speed_factor_pt  = to_double(t[i++]);
        pt_interval      = to_double(t[i++]);
        iterations       = (SQLINTEGER)to_long(t[i++]);
        iter_without_impr= (SQLINTEGER)to_long(t[i++]);
        max_removal_pct  = to_double(t[i++]);
        max_deterior_pct = to_double(t[i++]);
        clustered_removal= to_bit(t[i++]);
        idarp_adapted    = to_bit(t[i++]);
        without_operator = to_bit(t[i++]);
        terminals_per_node = (SQLINTEGER)to_long(t[i++]);

        priority_non_pt  = to_double(t[i++]);
        priority_pt      = to_double(t[i++]);
        priority_long_dist = to_double(t[i++]);
        priority_short_dist= to_double(t[i++]);

        current_iteration = (SQLINTEGER)to_long(t[i++]);
        total_distance    = to_double(t[i++]);
        computation_time  = to_double(t[i++]);

        avg_actual_ride_time      = to_double(t[i++]);
        avg_actual_ride_time_PT   = to_double(t[i++]);
        avg_actual_ride_time_noPT = to_double(t[i++]);

        avg_excess_ride_time      = to_double(t[i++]);
        avg_excess_ride_time_PT   = to_double(t[i++]);
        avg_excess_ride_time_noPT = to_double(t[i++]);

        avg_detour_factor         = to_double(t[i++]);
        avg_detour_factor_PT      = to_double(t[i++]);
        avg_detour_factor_noPT    = to_double(t[i++]);

        avg_waiting_time_PT       = to_double(t[i++]);
        pt_usage_ratio            = to_double(t[i++]);

        accept_random_removal     = (SQLINTEGER)to_long(t[i++]);
        accept_worst_removal      = (SQLINTEGER)to_long(t[i++]);
        accept_related_removal    = (SQLINTEGER)to_long(t[i++]);
        accept_route_removal      = (SQLINTEGER)to_long(t[i++]);

        accept_random_order_insertion = (SQLINTEGER)to_long(t[i++]);
        accept_greedy_insertion       = (SQLINTEGER)to_long(t[i++]);
        accept_two_regret_insertion   = (SQLINTEGER)to_long(t[i++]);

        improve_random_removal     = (SQLINTEGER)to_long(t[i++]);
        improve_worst_removal      = (SQLINTEGER)to_long(t[i++]);
        improve_related_removal    = (SQLINTEGER)to_long(t[i++]);
        improve_route_removal      = (SQLINTEGER)to_long(t[i++]);

        improve_random_order_insertion = (SQLINTEGER)to_long(t[i++]);
        improve_greedy_insertion       = (SQLINTEGER)to_long(t[i++]);
        improve_two_regret_insertion   = (SQLINTEGER)to_long(t[i++]);

        SQLRETURN rc = SQLExecute(stmt);
        if (rc != SQL_SUCCESS && rc != SQL_SUCCESS_WITH_INFO) {
            throw std::runtime_error("INSERT failed at line " + std::to_string(line_no));
        }

        // Important if you reuse the same stmt for multiple inserts (close cursor)
        SQLFreeStmt(stmt, SQL_CLOSE);
    };

            // Commit if everything worked
        check_rc(SQLExecDirectA(stmt, (SQLCHAR*)"COMMIT TRANSACTION", SQL_NTS), "COMMIT TRANSACTION");
    }

    catch (...) {
    SQLExecDirectA(stmt, (SQLCHAR*)"ROLLBACK TRANSACTION", SQL_NTS);
    throw;
}
    
};


void import_results_pipeline(){
    std::vector<std::string> data_darp = {"results_summary_datafile0.txt.txt", 
        "results_summary_datafile1.txt.txt", 
        "results_summary_datafile2.txt.txt", 
        "results_summary_datafile3.txt.txt",
        "results_summary_datafile4.txt.txt",
        "results_summary_datafile5.txt.txt",
        "results_summary_datafile6.txt.txt",
        "results_summary_datafile7.txt.txt",
        "results_summary_datafile8.txt.txt",
        "results_summary_datafile9.txt.txt",};

    std::string connStr = produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = connect_to_db(connStr, sql_param);

    create_results_table_if_missing(sql_param);

    sql_param = free_up_statement_for_new_query_type(sql_param);

    int instance_id = 1;

    for (const std::string& file : data_darp){
        std::string path = "C:\\Users\\enzot\\Documents\\Cesure\\1ere cesure inria Lille\\Codes\\Stage-Inria-01-09-2025--28-02-2026\\code_and_instances\\" + file ;
        import_results_summary_tsv(path, instance_id, sql_param);
        instance_id++ ;
    };
    sql_param = free_up_statement_for_new_query_type(sql_param);
    sql_param = cleanup_sql_query(sql_param);

};


// Interact with db for certain metrics
double avg_distance(int iteration, sql_params sql_param){
    double mean = 0;

    std::string sql = "SELECT AVG(total_distance) "
                        "FROM dbo.results_summary "
                        "WHERE current_iteration = ?;";

    SQLHSTMT stmt = sql_param.stmt;

    SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Bind instance_id parameter
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &iteration, 0, nullptr);

    // Execute query
    SQLExecute(stmt);

    // Bind result columns
    SQLBindCol(stmt, 1, SQL_C_DOUBLE,  &mean, 0, nullptr);

    if (SQLFetch(stmt) != SQL_SUCCESS) {
        throw std::runtime_error("No Mean for total distance found");
    };

    return mean;
};


double avg_ride_time(int iteration, sql_params sql_param){
    double mean = 0;

    std::string sql = "SELECT AVG(avg_actual_ride_time) "
                        "FROM dbo.results_summary "
                        "WHERE current_iteration = ?;";

    SQLHSTMT stmt = sql_param.stmt;

    SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Bind instance_id parameter
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &iteration, 0, nullptr);


    // Execute query
    SQLExecute(stmt);

    // Bind result columns
    SQLBindCol(stmt, 1, SQL_C_DOUBLE,  &mean, 0, nullptr);

    if (SQLFetch(stmt) != SQL_SUCCESS) {
        throw std::runtime_error("No Mean for ride time found");
    };

    return mean;
}


double avg_detour_factor(int iteration, sql_params sql_param){
    double mean = 0;

    std::string sql = "SELECT AVG(avg_detour_factor) "
                        "FROM dbo.results_summary "
                        "WHERE current_iteration = ?;";

    SQLHSTMT stmt = sql_param.stmt;

    SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Bind instance_id parameter
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &iteration, 0, nullptr);

    // Execute query
    SQLExecute(stmt);

    // Bind result columns
    SQLBindCol(stmt, 1, SQL_C_DOUBLE,  &mean, 0, nullptr);

    if (SQLFetch(stmt) != SQL_SUCCESS) {
        throw std::runtime_error("No Mean for detour factor found");
    };

    return mean;
}


void get_avg_pipeline(){
    int iteration = 25000;

    std::string connStr = produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = connect_to_db(connStr, sql_param);

    sql_param = free_up_statement_for_new_query_type(sql_param);

    double avg_dist = avg_distance(iteration, sql_param);

    sql_param = free_up_statement_for_new_query_type(sql_param);

    double avg_detour = avg_detour_factor(iteration, sql_param);

    sql_param = free_up_statement_for_new_query_type(sql_param);

    double avg_ride_time_val = avg_ride_time(iteration, sql_param);

    sql_param = free_up_statement_for_new_query_type(sql_param);

    sql_param = cleanup_sql_query(sql_param);

    std::cout << "The average distance travelled in all 10 instances is: " << avg_dist << "\n";
    std::cout << "The average detour factor in all 10 instances is: " << avg_detour << "\n";
    std::cout << "The average ride time in all 10 instances is: " << avg_ride_time_val << "\n";
}


// Main
int main(){
    get_avg_pipeline();
    return 1;
};










