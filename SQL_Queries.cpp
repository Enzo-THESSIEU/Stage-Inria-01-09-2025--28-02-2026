// Inclusions
#include "SQL_Queries.h"
#include "structures.h"

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


// Error Function
static void throw_odbc_error(const std::string& where,
                             SQLSMALLINT handleType,
                             SQLHANDLE handle)
{
    std::ostringstream oss;
    oss << "ODBC ERROR at " << where << "\n";

    SQLSMALLINT rec = 1;
    SQLCHAR sqlState[6] = {0};
    SQLCHAR msg[1024] = {0};
    SQLINTEGER nativeErr = 0;
    SQLSMALLINT msgLen = 0;

    bool hasDiag = false;

    while (SQLGetDiagRecA(handleType,
                          handle,
                          rec,
                          sqlState,
                          &nativeErr,
                          msg,
                          sizeof(msg),
                          &msgLen) == SQL_SUCCESS)
    {
        hasDiag = true;
        oss << "  [" << rec << "] "
            << "SQLSTATE=" << sqlState
            << ", NativeError=" << nativeErr
            << "\n      Message: " << msg << "\n";
        ++rec;
    }

    if (!hasDiag) {
        oss << "  <no diagnostic records available>\n";
    }

    throw std::runtime_error(oss.str());
}

void check_rc(SQLRETURN rc, const char* where) {
    if (rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO) return;
    throw std::runtime_error(std::string("ODBC error at: ") + where);
}


// TXT file reading functions
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


// Table Creating Functions
void table_creators::create_results_table_if_missing(const sql_params& sql_param)
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
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)){
        throw_odbc_error("SQLExecDirectA(create table)", SQL_HANDLE_STMT, stmt);
        }

};


void table_creators::create_request_table_if_missing(const sql_params& sql_param)
{
    SQLHSTMT stmt = sql_param.stmt;

    const std::string sql = R"SQL(
IF OBJECT_ID(N'dbo.requests', N'U') IS NULL
BEGIN
    CREATE TABLE dbo.requests (
        id INT IDENTITY(1,1) PRIMARY KEY,
        request_id INT, 
        instance_id INT NOT NULL,
        pickup_node_id INT,
        dropoff_node_id INT,

        t_pickup_dropoff FLOAT,
        mrt FLOAT,

        passengers INT,

        twL_pickup FLOAT,
        twU_pickup FLOAT,
        twL_dropoff FLOAT,
        twU_dropoff FLOAT,

        tw_score FLOAT,

        Correlated_request_1  FLOAT,
        Correlated_request_2  FLOAT,
        Correlated_request_3  FLOAT,
        Correlated_request_4  FLOAT,
        Correlated_request_5  FLOAT,
        Correlated_request_6  FLOAT,
        Correlated_request_7  FLOAT,
        Correlated_request_8  FLOAT,
        Correlated_request_9  FLOAT,
        Correlated_request_10 FLOAT,
        Correlated_request_11 FLOAT,
        Correlated_request_12 FLOAT,
        Correlated_request_13 FLOAT,
        Correlated_request_14 FLOAT,
        Correlated_request_15 FLOAT,
        Correlated_request_16 FLOAT,
        Correlated_request_17 FLOAT,
        Correlated_request_18 FLOAT,
        Correlated_request_19 FLOAT,
        Correlated_request_20 FLOAT,
        Correlated_request_21 FLOAT,
        Correlated_request_22 FLOAT,
        Correlated_request_23 FLOAT,
        Correlated_request_24 FLOAT,
        Correlated_request_25 FLOAT,
        Correlated_request_26 FLOAT,
        Correlated_request_27 FLOAT,
        Correlated_request_28 FLOAT,
        Correlated_request_29 FLOAT,
        Correlated_request_30 FLOAT,
        Correlated_request_31 FLOAT,
        Correlated_request_32 FLOAT,
        Correlated_request_33 FLOAT,
        Correlated_request_34 FLOAT,
        Correlated_request_35 FLOAT,
        Correlated_request_36 FLOAT,
        Correlated_request_37 FLOAT,
        Correlated_request_38 FLOAT,
        Correlated_request_39 FLOAT,
        Correlated_request_40 FLOAT,
        Correlated_request_41 FLOAT,
        Correlated_request_42 FLOAT,
        Correlated_request_43 FLOAT,
        Correlated_request_44 FLOAT,
        Correlated_request_45 FLOAT,
        Correlated_request_46 FLOAT,
        Correlated_request_47 FLOAT,
        Correlated_request_48 FLOAT,
        Correlated_request_49 FLOAT,
        Correlated_request_50 FLOAT,

        CONSTRAINT FK_correlated_requests_instance
        FOREIGN KEY (instance_id)
        REFERENCES dbo.instances(id)

    );
END
)SQL";

    SQLRETURN rc = SQLExecDirectA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)){
        throw_odbc_error("SQLExecDirectA(create table)", SQL_HANDLE_STMT, stmt);
        }

};


void table_creators::delete_table(const sql_params& sql_param, const std::string& table){
    // table should be trusted / not user input (identifiers cannot be parameterized)
    std::string sql =
        "IF OBJECT_ID(N'" + table + "', N'U') IS NOT NULL "
        "DROP TABLE " + table + ";";

    SQLRETURN rc = SQLExecDirectA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
        throw_odbc_error("SQLExecDirectA(delete_table)", SQL_HANDLE_STMT, sql_param.stmt);
    }
}

// Table Modifying Functions
void table_modifiers::add_column_to_table(const sql_params& sql_param,
                         const std::string& table,
                         const std::string& newColumn,
                         const std::string& columnType) {
    // NOTE: Identifiers (table/column/type) cannot be parameter-bound safely in ODBC.
    // Only call this with trusted strings you control.
    const std::string sql =
        "IF COL_LENGTH('" + table + "', '" + newColumn + "') IS NULL\n"
        "BEGIN\n"
        "    ALTER TABLE " + table + " ADD " + newColumn + " " + columnType + " NULL;\n"
        "END;";

    
    // Execute
    SQLRETURN rc = SQLExecDirectA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
        throw_odbc_error("SQLExecDirect", SQL_HANDLE_STMT, sql_param.stmt);
    }

    // Commit (ODBC default is autocommit ON; but commit explicitly is fine)
    SQLEndTran(SQL_HANDLE_DBC, sql_param.dbc, SQL_COMMIT);

    std::cout << "Column '" << newColumn << "' added to " << table
              << " (if it did not already exist)\n";
}



// Info Getting Functions (SELECT)

std::vector<int> info_getters::get_instances(sql_params sql_param) {
    std::vector<int> instance_ids;

    // IMPORTANT: keep spaces + add semicolon
    const std::string sql_req =
        "SELECT id "
        "FROM dbo.instances "
        "ORDER BY id;";

    SQLRETURN rc = SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_req.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLPrepare", SQL_HANDLE_STMT, sql_param.stmt);

    rc = SQLExecute(sql_param.stmt);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLExecute", SQL_HANDLE_STMT, sql_param.stmt);

    // 4) Bind result column (id) + fetch all rows
    SQLINTEGER id = 0;
    SQLLEN id_ind = 0;

    rc = SQLBindCol(sql_param.stmt, 1, SQL_C_SLONG, &id, 0, &id_ind);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindCol(id)", SQL_HANDLE_STMT, sql_param.stmt);

    while (true) {
        rc = SQLFetch(sql_param.stmt);
        if (rc == SQL_NO_DATA) break; // normal end of result set
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLFetch", SQL_HANDLE_STMT, sql_param.stmt);

        if (id_ind != SQL_NULL_DATA)
            instance_ids.push_back((int)id);
    }

    // Optional: treat empty table as error (remove if you prefer empty vector)
    if (instance_ids.empty())
        throw std::runtime_error("No rows returned");

    return instance_ids;
};



instance_info info_getters::get_instance_info(int instance_id, sql_params sql_param){
    instance_info instance{};
    instance.instance_id = instance_id;

    std::string sql_req = 
    "SELECT n_vehicles, n_requests, n_terminals "
    "FROM dbo.instances "
    "WHERE id = ?";

    SQLRETURN rc = SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_req.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLPrepare (get_instance_info) ", SQL_HANDLE_STMT, sql_param.stmt);

    // Bind instance_id parameter
    SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &instance_id, 0, nullptr);

    // Execute query
    SQLExecute(sql_param.stmt);

    // Bind result columns
    SQLBindCol(sql_param.stmt, 1, SQL_C_SLONG,  &instance.n_vehicles, 0, nullptr);
    SQLBindCol(sql_param.stmt, 2, SQL_C_SLONG,  &instance.n_requests,       0, nullptr);
    SQLBindCol(sql_param.stmt, 3, SQL_C_SLONG,  &instance.n_terminals,       0, nullptr);

    // Fetch first pickup row
    rc = SQLFetch(sql_param.stmt);
        if (rc == SQL_NO_DATA) {
            throw std::runtime_error(
                "No information found for instance " + std::to_string(instance_id)
            );
        }
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
            throw_odbc_error(
                "SQLFetch(get_instance_info) failed for instance " + std::to_string(instance_id),
                SQL_HANDLE_STMT, sql_param.stmt
            );
        };

    return instance;
};



std::pair<std::vector<node_info>, std::vector<node_info>> info_getters::link_pickups_dropoffs(int instance_id, instance_info instance, sql_params sql_param)
{
    std::vector<node_info> pickup_nodes = {};
    std::vector<node_info> dropoff_nodes = {};

    /* ============================================================
       1. SQL queries
       ============================================================ */

    const std::string sql_pickup =
        "SELECT node_id, x, y, twL, twU "
        "FROM dbo.nodes "
        "WHERE instance_id = ? AND node_id = ? AND node_type = 'Pickup Node';";

    const std::string sql_dropoff =
        "SELECT node_id, x, y, twL, twU "
        "FROM dbo.nodes "
        "WHERE instance_id = ? AND node_id = ? AND node_type = 'Dropoff Node';";

    /* ============================================================
       3. Fetch PICKUP node
       ============================================================ */

    for (int pickup_node_id = 0; pickup_node_id < instance.n_requests; ++pickup_node_id){

        node_info pickup{};

        SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_pickup.c_str(), SQL_NTS);

        // Bind instance_id parameter
        SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &instance_id, 0, nullptr);

        SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &pickup_node_id, 0, nullptr);

        // Execute query
        SQLExecute(sql_param.stmt);

        // Bind result columns
        SQLBindCol(sql_param.stmt, 1, SQL_C_SLONG,  &pickup.node_id, 0, nullptr);
        SQLBindCol(sql_param.stmt, 2, SQL_C_DOUBLE, &pickup.x,       0, nullptr);
        SQLBindCol(sql_param.stmt, 3, SQL_C_DOUBLE, &pickup.y,       0, nullptr);
        SQLBindCol(sql_param.stmt, 4, SQL_C_DOUBLE, &pickup.twL,     0, nullptr);
        SQLBindCol(sql_param.stmt, 5, SQL_C_DOUBLE, &pickup.twU,     0, nullptr);

        // Fetch first pickup row
        if (SQLFetch(sql_param.stmt) != SQL_SUCCESS) {
            throw std::runtime_error("No pickup node found");
        }

        SQLFreeStmt(sql_param.stmt, SQL_CLOSE);
        SQLFreeStmt(sql_param.stmt, SQL_UNBIND);
        SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS);


    /* ============================================================
       4. Fetch DROPOFF node corresponding to pickup
       ============================================================ */

        node_info dropoff{};
        int dropoff_node_id = pickup_node_id + instance.n_requests;

        SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_dropoff.c_str(), SQL_NTS);

        // Bind parameters
        SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &instance_id, 0, nullptr);

        SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &dropoff_node_id, 0, nullptr);

        // Execute query
        SQLExecute(sql_param.stmt);

        // Bind result columns
        SQLBindCol(sql_param.stmt, 1, SQL_C_SLONG,  &dropoff.node_id, 0, nullptr);
        SQLBindCol(sql_param.stmt, 2, SQL_C_DOUBLE, &dropoff.x,       0, nullptr);
        SQLBindCol(sql_param.stmt, 3, SQL_C_DOUBLE, &dropoff.y,       0, nullptr);
        SQLBindCol(sql_param.stmt, 4, SQL_C_DOUBLE, &dropoff.twL,     0, nullptr);
        SQLBindCol(sql_param.stmt, 5, SQL_C_DOUBLE, &dropoff.twU,     0, nullptr);

        if (SQLFetch(sql_param.stmt) != SQL_SUCCESS) {
            throw std::runtime_error("No dropoff node found");
        }

        SQLFreeStmt(sql_param.stmt, SQL_CLOSE);
        SQLFreeStmt(sql_param.stmt, SQL_UNBIND);
        SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS);


        pickup_nodes.push_back(pickup);
        dropoff_nodes.push_back(dropoff);

    };

    /* ============================================================
       6. Return both nodes
       ============================================================ */

    return { pickup_nodes, dropoff_nodes };
}



node_info info_getters::get_node_info(int instance_id, int node_id, sql_params sql_param){
    node_info node{};

    SQLHSTMT stmt = sql_param.stmt;

    std::string sql = "SELECT node_id, x, y, twL, twU, mobile, wheelchair, service_dur, node_type, twL_constrained, twU_constrained "
                        "FROM dbo.nodes "
                        "WHERE instance_id = ? AND node_id =?;";

    SQLRETURN rc = SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLPrepareA(get_node_info)", SQL_HANDLE_STMT, stmt);

        rc = SQLBindParameter(stmt, 1, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLBindParameter(instance_id)", SQL_HANDLE_STMT, stmt);

        rc = SQLBindParameter(stmt, 2, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &node_id, 0, nullptr);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLBindParameter(node_id)", SQL_HANDLE_STMT, stmt);

    // Bind Results

    char node_type_buf[32];   // Necessary since odbc can't bind to string directly

    SQLBindCol(stmt, 1, SQL_C_SLONG,  &node.node_id,        0, nullptr);
    SQLBindCol(stmt, 2, SQL_C_DOUBLE, &node.x,              0, nullptr);
    SQLBindCol(stmt, 3, SQL_C_DOUBLE, &node.y,              0, nullptr);
    SQLBindCol(stmt, 4, SQL_C_DOUBLE, &node.twL,            0, nullptr);
    SQLBindCol(stmt, 5, SQL_C_DOUBLE, &node.twU,            0, nullptr);
    SQLBindCol(stmt, 6, SQL_C_SLONG,  &node.mobile,         0, nullptr);
    SQLBindCol(stmt, 7, SQL_C_SLONG,  &node.wheelchair,     0, nullptr);
    SQLBindCol(stmt, 8, SQL_C_DOUBLE, &node.service_dur,    0, nullptr);
    SQLBindCol(stmt, 9, SQL_C_CHAR,   node_type_buf,        sizeof(node_type_buf), nullptr);
    SQLBindCol(stmt,10, SQL_C_DOUBLE, &node.twL_constrained,0, nullptr);
    SQLBindCol(stmt,11, SQL_C_DOUBLE, &node.twU_constrained,0, nullptr);

    rc = SQLExecute(stmt);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLExecute(get_node_info)", SQL_HANDLE_STMT, stmt);

    rc = SQLFetch(stmt);
    if (rc == SQL_NO_DATA) {
        throw std::runtime_error(
            "Node " + std::to_string(node_id) + " not found for instance " + std::to_string(instance_id)
        );
    }
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
        throw_odbc_error(
            "SQLFetch(get_node_info) failed for node " + std::to_string(node_id) +
            " instance " + std::to_string(instance_id),
            SQL_HANDLE_STMT, stmt
        );
    };

    node.node_type = node_type_buf;

    return node;


};

std::vector<db_request_info> info_getters::get_n_requests(int n, int offset, int instance_id, sql_params sql_param)
{
    const std::string sql = R"SQL(
        SELECT
            request_id,
            instance_id,
            pickup_node_id,
            dropoff_node_id,
            t_pickup_dropoff,
            mrt,
            passengers,
            twL_pickup,
            twU_pickup,
            twL_dropoff,
            twU_dropoff
        FROM dbo.requests
        WHERE instance_id = ?
        AND tw_score IS NOT NULL
        ORDER BY tw_score ASC
        OFFSET ? ROWS
        FETCH NEXT ? ROWS ONLY;
    )SQL";

    SQLHSTMT stmt = sql_param.stmt;

    SQLRETURN rc = SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLPrepareA(get_n_requests)", SQL_HANDLE_STMT, stmt);

    // ---- INPUT params ----
    rc = SQLBindParameter(stmt, 3, SQL_PARAM_INPUT,
                          SQL_C_SLONG, SQL_INTEGER, 0, 0,
                          &n, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(n)", SQL_HANDLE_STMT, stmt);

    rc = SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                          SQL_C_SLONG, SQL_INTEGER, 0, 0,
                          &instance_id, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(instance_id)", SQL_HANDLE_STMT, stmt);

    rc = SQLBindParameter(stmt, 2, SQL_PARAM_INPUT,
                          SQL_C_SLONG, SQL_INTEGER, 0, 0,
                          &offset, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(offset)", SQL_HANDLE_STMT, stmt);

    // ---- OUTPUT cols ----
    db_request_info req{};

    // Indicators (IMPORTANT)
    SQLLEN ind_request_id=0, ind_instance_id=0, ind_pickup=0, ind_dropoff=0;
    SQLLEN ind_tpd=0, ind_mrt=0, ind_passengers=0;
    SQLLEN ind_twLp=0, ind_twUp=0, ind_twLd=0, ind_twUd=0;

    auto bindcol = [&](SQLUSMALLINT col, SQLSMALLINT c_type, void* buf, SQLLEN buflen, SQLLEN* ind, const char* where){
        SQLRETURN r = SQLBindCol(stmt, col, c_type, buf, buflen, ind);
        if (!(r == SQL_SUCCESS || r == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error(where, SQL_HANDLE_STMT, stmt);
    };

    bindcol(1,  SQL_C_SLONG,  &req.request_id,       0, &ind_request_id, "SQLBindCol(request_id)");
    bindcol(2,  SQL_C_SLONG,  &req.instance_id,      0, &ind_instance_id,"SQLBindCol(instance_id)");
    bindcol(3,  SQL_C_SLONG,  &req.pickup_node_id,   0, &ind_pickup,     "SQLBindCol(pickup_node_id)");
    bindcol(4,  SQL_C_SLONG,  &req.dropoff_node_id,  0, &ind_dropoff,    "SQLBindCol(dropoff_node_id)");
    bindcol(5,  SQL_C_DOUBLE, &req.t_pickup_dropoff, 0, &ind_tpd,        "SQLBindCol(t_pickup_dropoff)");
    bindcol(6,  SQL_C_DOUBLE, &req.mrt,              0, &ind_mrt,        "SQLBindCol(mrt)");
    bindcol(7,  SQL_C_SLONG,  &req.passengers,       0, &ind_passengers, "SQLBindCol(passengers)");
    bindcol(8,  SQL_C_DOUBLE, &req.twL_pickup,       0, &ind_twLp,       "SQLBindCol(twL_pickup)");
    bindcol(9,  SQL_C_DOUBLE, &req.twU_pickup,       0, &ind_twUp,       "SQLBindCol(twU_pickup)");
    bindcol(10, SQL_C_DOUBLE, &req.twL_dropoff,      0, &ind_twLd,       "SQLBindCol(twL_dropoff)");
    bindcol(11, SQL_C_DOUBLE, &req.twU_dropoff,      0, &ind_twUd,       "SQLBindCol(twU_dropoff)");

    // ---- Execute + fetch ----
    rc = SQLExecute(stmt);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLExecute(get_n_requests)", SQL_HANDLE_STMT, stmt);

    std::vector<db_request_info> out;
    out.reserve(n);

    while ((rc = SQLFetch(stmt)) != SQL_NO_DATA) {
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLFetch(get_n_requests)", SQL_HANDLE_STMT, stmt);

        // Optional: handle NULLs safely (set defaults)
        if (ind_pickup   == SQL_NULL_DATA) req.pickup_node_id  = -1;
        if (ind_dropoff  == SQL_NULL_DATA) req.dropoff_node_id = -1;
        if (ind_tpd      == SQL_NULL_DATA) req.t_pickup_dropoff = 0.0;
        if (ind_mrt      == SQL_NULL_DATA) req.mrt = 0.0;
        if (ind_twLp     == SQL_NULL_DATA) req.twL_pickup = 0.0;
        if (ind_twUp     == SQL_NULL_DATA) req.twU_pickup = 0.0;
        if (ind_twLd     == SQL_NULL_DATA) req.twL_dropoff = 0.0;
        if (ind_twUd     == SQL_NULL_DATA) req.twU_dropoff = 0.0;

        out.push_back(req);
    }

    return out;
}

std::vector<node_info> info_getters::get_nodes_by_type(int instance_id, const std::string& node_type, int limit, sql_params sql_param)
{
    std::vector<node_info> out;
    out.reserve(limit);

    SQLHSTMT stmt = sql_param.stmt;

    const std::string sql = R"SQL(
        SELECT TOP (?)
            node_id, x, y, twL, twU, mobile, wheelchair, service_dur, node_type, twL_constrained, twU_constrained
        FROM dbo.nodes
        WHERE instance_id = ?
          AND node_type = ?
        ORDER BY node_id ASC;
    )SQL";

    SQLRETURN rc = SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLPrepareA(get_nodes_by_type)", SQL_HANDLE_STMT, stmt);

    // inputs
    rc = SQLBindParameter(stmt, 1, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &limit, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(limit)", SQL_HANDLE_STMT, stmt);

    rc = SQLBindParameter(stmt, 2, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(instance_id)", SQL_HANDLE_STMT, stmt);

    // node_type buffer
    char node_type_in[64]{};
    std::snprintf(node_type_in, sizeof(node_type_in), "%s", node_type.c_str());
    rc = SQLBindParameter(stmt, 3, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_VARCHAR, 0, 0, node_type_in, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(node_type)", SQL_HANDLE_STMT, stmt);

    // outputs
    node_info node{};
    char node_type_buf[64]{};

    SQLLEN ind[11]{};

    auto bindcol = [&](SQLUSMALLINT col, SQLSMALLINT ctype, void* buf, SQLLEN buflen, SQLLEN* indp, const char* where){
        SQLRETURN r = SQLBindCol(stmt, col, ctype, buf, buflen, indp);
        if (!(r == SQL_SUCCESS || r == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error(where, SQL_HANDLE_STMT, stmt);
    };

    bindcol(1,  SQL_C_SLONG,  &node.node_id,         0, &ind[0],  "SQLBindCol(node_id)");
    bindcol(2,  SQL_C_DOUBLE, &node.x,               0, &ind[1],  "SQLBindCol(x)");
    bindcol(3,  SQL_C_DOUBLE, &node.y,               0, &ind[2],  "SQLBindCol(y)");
    bindcol(4,  SQL_C_DOUBLE, &node.twL,             0, &ind[3],  "SQLBindCol(twL)");
    bindcol(5,  SQL_C_DOUBLE, &node.twU,             0, &ind[4],  "SQLBindCol(twU)");
    bindcol(6,  SQL_C_SLONG,  &node.mobile,          0, &ind[5],  "SQLBindCol(mobile)");
    bindcol(7,  SQL_C_SLONG,  &node.wheelchair,      0, &ind[6],  "SQLBindCol(wheelchair)");
    bindcol(8,  SQL_C_DOUBLE, &node.service_dur,     0, &ind[7],  "SQLBindCol(service_dur)");
    bindcol(9,  SQL_C_CHAR,   node_type_buf, sizeof(node_type_buf), &ind[8],  "SQLBindCol(node_type)");
    bindcol(10, SQL_C_DOUBLE, &node.twL_constrained, 0, &ind[9],  "SQLBindCol(twL_constrained)");
    bindcol(11, SQL_C_DOUBLE, &node.twU_constrained, 0, &ind[10], "SQLBindCol(twU_constrained)");

    rc = SQLExecute(stmt);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLExecute(get_nodes_by_type)", SQL_HANDLE_STMT, stmt);

    while ((rc = SQLFetch(stmt)) != SQL_NO_DATA) {
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLFetch(get_nodes_by_type)", SQL_HANDLE_STMT, stmt);

        node.node_type = node_type_buf;
        out.push_back(node);

        if ((int)out.size() >= limit) break;
    }

    return out;
}



double info_getters::avg_distance(int iteration, sql_params sql_param){
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

double info_getters::avg_ride_time(int iteration, sql_params sql_param){
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

double info_getters::avg_detour_factor(int iteration, sql_params sql_param){
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



// Input to DB Functions
void input_to_db::import_results_summary_tsv(const std::string& path, int instance_id, sql_params& sql_param)
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



void input_to_db::input_constrained_TW_into_db(std::vector<node_info> pickup_nodes, std::vector<node_info> dropoff_nodes, int instance_id, sql_params& sql_param){
    std::string sql_update_pickup = 
    "UPDATE dbo.nodes "
    "SET twL_constrained = ?, twU_constrained = ? "
    "WHERE instance_id = ? AND node_id = ?";

    std::string sql_update_dropoff = 
    "UPDATE dbo.nodes "
    "SET twL_constrained = ?, twU_constrained = ? "
    "WHERE instance_id = ? AND node_id = ?";

    /* ============================================================
    5. Interact with the DB
    ============================================================ */

    // Update Pickup Windows

    for (node_info pickup_node : pickup_nodes){

    SQLAllocHandle(SQL_HANDLE_STMT, sql_param.dbc, &sql_param.stmt);
    SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_update_pickup.c_str(), SQL_NTS);

    // Bind parametres
    SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &pickup_node.twL_constrained, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &pickup_node.twU_constrained, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 3, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 4, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &pickup_node.node_id, 0, nullptr);

    // Execute Query
    SQLExecute(sql_param.stmt);

    }

    SQLFreeStmt(sql_param.stmt, SQL_CLOSE);
    SQLFreeStmt(sql_param.stmt, SQL_UNBIND);
    SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS);

    for (node_info dropoff_node : dropoff_nodes){

    // Update Dropoff Windows
    SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql_update_dropoff.c_str(), SQL_NTS);

    // Bind parametres
    SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &dropoff_node.twL_constrained, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &dropoff_node.twU_constrained, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 3, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 4, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &dropoff_node.node_id, 0, nullptr);

    // Execute Query
    SQLExecute(sql_param.stmt);

};

    SQLFreeStmt(sql_param.stmt, SQL_CLOSE);
    SQLFreeStmt(sql_param.stmt, SQL_UNBIND);
    SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS);

};



int input_to_db::insert_trains_into_db(const std::vector<Train>& all_trains, sql_params& sql_param) {

    // Prepare insert
    std::string sql =
        "INSERT INTO dbo.Timetable (TripName, DepartureTime, ArrivalTime) "
        "VALUES (?, ?, ?);";

    SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Buffers for parameters (reused each loop)
    char stations[256];
    double departure;
    double arrival;

    // Bind once (buffers updated per row)
    SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_VARCHAR, sizeof(stations) - 1, 0, stations, sizeof(stations), nullptr);
    SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &departure, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 3, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &arrival,  0, nullptr);

    for (const auto& train : all_trains){
        // Times
        std::string stations = train.stations;
        int departure = train.departure;
        int arrival = train.arrival;


        SQLRETURN rc = SQLExecute(sql_param.stmt);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
            // rollback transaction properly
            SQLEndTran(SQL_HANDLE_DBC, sql_param.dbc, SQL_ROLLBACK);
            throw_odbc_error("SQLExecute(insert)", SQL_HANDLE_STMT, sql_param.stmt);
        }

        SQLFreeStmt(sql_param.stmt, SQL_CLOSE);
    }


    // Commit at the end (only matters if autocommit is OFF)
    SQLEndTran(SQL_HANDLE_DBC, sql_param.dbc, SQL_COMMIT);

    // Optional: also reset bindings/params if this stmt will be reused for other queries
    SQLFreeStmt(sql_param.stmt, SQL_UNBIND);
    SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS);

    return 0;
}



void input_to_db::insert_correlated_set(const sql_params& sql_param, std::vector<int> correlated_set, int instance_id, int request_id){
std::string sql = R"SQL(
UPDATE dbo.requests
SET
    Correlated_request_1  = ?,  Correlated_request_2  = ?,  Correlated_request_3  = ?,  Correlated_request_4  = ?,  Correlated_request_5  = ?,
    Correlated_request_6  = ?,  Correlated_request_7  = ?,  Correlated_request_8  = ?,  Correlated_request_9  = ?,  Correlated_request_10 = ?,
    Correlated_request_11 = ?,  Correlated_request_12 = ?,  Correlated_request_13 = ?,  Correlated_request_14 = ?,  Correlated_request_15 = ?,
    Correlated_request_16 = ?,  Correlated_request_17 = ?,  Correlated_request_18 = ?,  Correlated_request_19 = ?,  Correlated_request_20 = ?,
    Correlated_request_21 = ?,  Correlated_request_22 = ?,  Correlated_request_23 = ?,  Correlated_request_24 = ?,  Correlated_request_25 = ?,
    Correlated_request_26 = ?,  Correlated_request_27 = ?,  Correlated_request_28 = ?,  Correlated_request_29 = ?,  Correlated_request_30 = ?,
    Correlated_request_31 = ?,  Correlated_request_32 = ?,  Correlated_request_33 = ?,  Correlated_request_34 = ?,  Correlated_request_35 = ?,
    Correlated_request_36 = ?,  Correlated_request_37 = ?,  Correlated_request_38 = ?,  Correlated_request_39 = ?,  Correlated_request_40 = ?,
    Correlated_request_41 = ?,  Correlated_request_42 = ?,  Correlated_request_43 = ?,  Correlated_request_44 = ?,  Correlated_request_45 = ?,
    Correlated_request_46 = ?,  Correlated_request_47 = ?,  Correlated_request_48 = ?,  Correlated_request_49 = ?,  Correlated_request_50 = ?
WHERE instance_id = ? AND request_id = ?;
)SQL";

    check_rc(SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS), "SQLPrepareA");

    for (int i = 1; i < 51; ++i){
        int value = correlated_set[i];
        SQLBindParameter(sql_param.stmt, i, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &value, 0, nullptr);
    };

    SQLBindParameter(sql_param.stmt, 51, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
    SQLBindParameter(sql_param.stmt, 52, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &request_id, 0, nullptr);

    SQLRETURN rc = SQLExecute(sql_param.stmt);
        if (rc != SQL_SUCCESS && rc != SQL_SUCCESS_WITH_INFO) {
            throw std::runtime_error("UPDATE dbo.requests failed");
        };
}



void input_to_db::input_to_db_tw_scores(const sql_params& sql_param, double tw_score, int instance_id, int request_id){
    std::string sql = R"SQL(
    UPDATE dbo.requests
    SET tw_score = ?
    WHERE instance_id = ? AND request_id = ?;
    )SQL";


    

    SQLRETURN rc = SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLPrepareA(input_to_db_tw_scores)", SQL_HANDLE_STMT, sql_param.stmt);

    rc = SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT,
                          SQL_C_SLONG, SQL_INTEGER, 0, 0,
                          &instance_id, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(instance_id)", SQL_HANDLE_STMT, sql_param.stmt);

    rc = SQLBindParameter(sql_param.stmt, 3, SQL_PARAM_INPUT,
                          SQL_C_SLONG, SQL_INTEGER, 0, 0,
                          &request_id, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(request_id)", SQL_HANDLE_STMT, sql_param.stmt);

    rc = SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT,
                          SQL_C_DOUBLE, SQL_DOUBLE, 0, 0,
                          &tw_score, 0, nullptr);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLBindParameter(tw_score)", SQL_HANDLE_STMT, sql_param.stmt);

    rc = SQLExecute(sql_param.stmt);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLExecute(input_to_db_tw_scores)", SQL_HANDLE_STMT, sql_param.stmt);

    SQLLEN rows = 0;
    rc = SQLRowCount(sql_param.stmt, &rows);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
        throw_odbc_error("SQLRowCount(input_to_db_tw_scores)", SQL_HANDLE_STMT, sql_param.stmt);

    if (rows == 0) {
        throw std::runtime_error(
            "UPDATE updated 0 rows for request " + std::to_string(request_id) +
            " instance " + std::to_string(instance_id) + " (check keys exist)"
        );
    }

};



void input_to_db::input_to_db_request_info(const sql_params& sql_param, db_request_info req){

    std::string sql = R"SQL(
    INSERT INTO dbo.requests (
        request_id,
        instance_id,
        pickup_node_id,
        dropoff_node_id,
        t_pickup_dropoff,
        mrt,
        passengers,
        twL_pickup,
        twU_pickup,
        twL_dropoff,
        twU_dropoff
    )
    VALUES (
        ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
    );
    )SQL";


SQLRETURN rc = SQLPrepareA(sql_param.stmt, (SQLCHAR*)sql.c_str(), SQL_NTS); 
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) 
    throw_odbc_error("SQLPrepareA(input_to_db_tw_scores)", SQL_HANDLE_STMT, sql_param.stmt);

// 1) request_id
rc = SQLBindParameter(sql_param.stmt, 1, SQL_PARAM_INPUT,
                      SQL_C_SLONG, SQL_INTEGER, 0, 0,
                      &req.request_id, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(request_id)", SQL_HANDLE_STMT, sql_param.stmt);

// 2) instance_id
rc = SQLBindParameter(sql_param.stmt, 2, SQL_PARAM_INPUT,
                      SQL_C_SLONG, SQL_INTEGER, 0, 0,
                      &req.instance_id, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(instance_id)", SQL_HANDLE_STMT, sql_param.stmt);

// 3) pickup_node_id
rc = SQLBindParameter(sql_param.stmt, 3, SQL_PARAM_INPUT,
                      SQL_C_SLONG, SQL_INTEGER, 0, 0,
                      &req.pickup_node_id, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(pickup_node_id)", SQL_HANDLE_STMT, sql_param.stmt);

// 4) dropoff_node_id
rc = SQLBindParameter(sql_param.stmt, 4, SQL_PARAM_INPUT,
                      SQL_C_SLONG, SQL_INTEGER, 0, 0,
                      &req.dropoff_node_id, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(dropoff_node_id)", SQL_HANDLE_STMT, sql_param.stmt);

// 5) t_pickup_dropoff
rc = SQLBindParameter(sql_param.stmt, 5, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.t_pickup_dropoff, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(t_pickup_dropoff)", SQL_HANDLE_STMT, sql_param.stmt);

// 6) mrt
rc = SQLBindParameter(sql_param.stmt, 6, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.mrt, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(mrt)", SQL_HANDLE_STMT, sql_param.stmt);

// 7) passengers
rc = SQLBindParameter(sql_param.stmt, 7, SQL_PARAM_INPUT,
                      SQL_C_SLONG, SQL_INTEGER, 0, 0,
                      &req.passengers, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(passengers)", SQL_HANDLE_STMT, sql_param.stmt);

// 8) twL_pickup
rc = SQLBindParameter(sql_param.stmt, 8, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.twL_pickup, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(twL_pickup)", SQL_HANDLE_STMT, sql_param.stmt);

// 9) twU_pickup
rc = SQLBindParameter(sql_param.stmt, 9, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.twU_pickup, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(twU_pickup)", SQL_HANDLE_STMT, sql_param.stmt);

// 10) twL_dropoff
rc = SQLBindParameter(sql_param.stmt, 10, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.twL_dropoff, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(twL_dropoff)", SQL_HANDLE_STMT, sql_param.stmt);

// 11) twU_dropoff
rc = SQLBindParameter(sql_param.stmt, 11, SQL_PARAM_INPUT,
                      SQL_C_DOUBLE, SQL_FLOAT, 0, 0,
                      &req.twU_dropoff, 0, nullptr);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLBindParameter(twU_dropoff)", SQL_HANDLE_STMT, sql_param.stmt);



    rc = SQLExecute(sql_param.stmt);
if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
    throw_odbc_error("SQLExecute(insert dbo.requests)", SQL_HANDLE_STMT, sql_param.stmt);

};




