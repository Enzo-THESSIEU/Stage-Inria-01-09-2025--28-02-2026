#include <windows.h>
#include <sqlext.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

// g++ -g "Complete_time_windows.cpp" -o "Complete_time_windows.exe" -lodbc32 -lodbccp32


struct node_info{
        int node_id;
        double x;
        double y;
        double twL;
        double twU;
    };

struct instance_info{
    int n_vehicles;
    int n_requests;
    int n_terminals;
};


std::string produce_conn_str(const std::string& server, const std::string& database) {
    // Note: Encrypt=no for local dev; if you want Encrypt=yes, add TrustServerCertificate=yes
    return "Driver={ODBC Driver 18 for SQL Server};"
           "Server=" + server + ";"
           "Database=" + database + ";"
           "Trusted_Connection=yes;"
           "Encrypt=no;";
}


static void throw_odbc_error(const std::string& where, SQLSMALLINT handleType, SQLHANDLE handle) {
    SQLCHAR sqlState[6]{}, msg[1024]{};
    SQLINTEGER nativeErr = 0;
    SQLSMALLINT msgLen = 0;

    if (SQLGetDiagRecA(handleType, handle, 1, sqlState, &nativeErr, msg, sizeof(msg), &msgLen) == SQL_SUCCESS) {
        throw std::runtime_error(where + " | " + (char*)sqlState + " | " + (char*)msg);
    }
    throw std::runtime_error(where + " | ODBC error (no diagnostics)");
}


void add_column_to_table(const std::string& connStr,
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

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    // Allocate environment
    if (SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env) != SQL_SUCCESS)
        throw std::runtime_error("SQLAllocHandle ENV failed");

    if (SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0) != SQL_SUCCESS) {
        SQLFreeHandle(SQL_HANDLE_ENV, env);
        throw std::runtime_error("SQLSetEnvAttr ODBC version failed");
    }

    // Allocate connection
    if (SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc) != SQL_SUCCESS) {
        SQLFreeHandle(SQL_HANDLE_ENV, env);
        throw std::runtime_error("SQLAllocHandle DBC failed");
    }

    // Connect using connection string
    SQLCHAR outConnStr[2048];
    SQLSMALLINT outLen = 0;
    SQLRETURN rc = SQLDriverConnectA(
        dbc,
        nullptr,
        (SQLCHAR*)connStr.c_str(),
        SQL_NTS,
        outConnStr,
        sizeof(outConnStr),
        &outLen,
        SQL_DRIVER_NOPROMPT
    );
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
        throw_odbc_error("SQLDriverConnect", SQL_HANDLE_DBC, dbc);
    }

    // Allocate statement
    if (SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt) != SQL_SUCCESS) {
        SQLDisconnect(dbc);
        SQLFreeHandle(SQL_HANDLE_DBC, dbc);
        SQLFreeHandle(SQL_HANDLE_ENV, env);
        throw std::runtime_error("SQLAllocHandle STMT failed");
    }

    // Execute
    rc = SQLExecDirectA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);
    if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO)) {
        throw_odbc_error("SQLExecDirect", SQL_HANDLE_STMT, stmt);
    }

    // Commit (ODBC default is autocommit ON; but commit explicitly is fine)
    SQLEndTran(SQL_HANDLE_DBC, dbc, SQL_COMMIT);

    // Cleanup
    SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);

    std::cout << "Column '" << newColumn << "' added to " << table
              << " (if it did not already exist)\n";
}


std::vector<int> get_instances(const std::string& connStr) {
    std::vector<int> instance_ids;

    // IMPORTANT: keep spaces + add semicolon
    const std::string sql_req =
        "SELECT id "
        "FROM dbo.instances "
        "ORDER BY id;";

    SQLHENV env  = SQL_NULL_HENV;
    SQLHDBC dbc  = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    try {
        // 1) Allocate environment + set ODBC version
        SQLRETURN rc = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw std::runtime_error("SQLAllocHandle(ENV) failed");

        rc = SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw std::runtime_error("SQLSetEnvAttr(ODBC3) failed");

        // 2) Allocate connection + connect
        rc = SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw std::runtime_error("SQLAllocHandle(DBC) failed");

        rc = SQLDriverConnectA(
            dbc, nullptr,
            (SQLCHAR*)connStr.c_str(), SQL_NTS,
            nullptr, 0, nullptr,
            SQL_DRIVER_NOPROMPT
        );
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLDriverConnect", SQL_HANDLE_DBC, dbc);

        // 3) Prepare + execute query
        rc = SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw std::runtime_error("SQLAllocHandle(STMT) failed");

        rc = SQLPrepareA(stmt, (SQLCHAR*)sql_req.c_str(), SQL_NTS);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLPrepare", SQL_HANDLE_STMT, stmt);

        rc = SQLExecute(stmt);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLExecute", SQL_HANDLE_STMT, stmt);

        // 4) Bind result column (id) + fetch all rows
        SQLINTEGER id = 0;
        SQLLEN id_ind = 0;

        rc = SQLBindCol(stmt, 1, SQL_C_SLONG, &id, 0, &id_ind);
        if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
            throw_odbc_error("SQLBindCol(id)", SQL_HANDLE_STMT, stmt);

        while (true) {
            rc = SQLFetch(stmt);
            if (rc == SQL_NO_DATA) break; // normal end of result set
            if (!(rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO))
                throw_odbc_error("SQLFetch", SQL_HANDLE_STMT, stmt);

            if (id_ind != SQL_NULL_DATA)
                instance_ids.push_back((int)id);
        }

        // Optional: treat empty table as error (remove if you prefer empty vector)
        if (instance_ids.empty())
            throw std::runtime_error("No rows returned");

        // 5) Cleanup
        SQLFreeHandle(SQL_HANDLE_STMT, stmt);
        SQLDisconnect(dbc);
        SQLFreeHandle(SQL_HANDLE_DBC, dbc);
        SQLFreeHandle(SQL_HANDLE_ENV, env);

        return instance_ids;
    }
    catch (...) {
        // Cleanup on error too
        if (stmt != SQL_NULL_HSTMT) SQLFreeHandle(SQL_HANDLE_STMT, stmt);
        if (dbc  != SQL_NULL_HDBC)  { SQLDisconnect(dbc); SQLFreeHandle(SQL_HANDLE_DBC, dbc); }
        if (env  != SQL_NULL_HENV)  SQLFreeHandle(SQL_HANDLE_ENV, env);
        throw;
    }
}


instance_info get_instance_node_info(std::string connStr, int instance_id){
    instance_info instance{};

    std::string sql_req = 
    "SELECT n_vehicles, n_requests, n_terminals "
    "FROM dbo.instances "
    "WHERE id = ?";

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    // Allocate environment
    SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
    SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);

    // Allocate connection
    SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);

    // Connect using connection string
    SQLDriverConnectA(
        dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );

    SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
    SQLPrepareA(stmt, (SQLCHAR*)sql_req.c_str(), SQL_NTS);

    // Bind instance_id parameter
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &instance_id, 0, nullptr);

    // Execute query
    SQLExecute(stmt);

    // Bind result columns
    SQLBindCol(stmt, 1, SQL_C_SLONG,  &instance.n_vehicles, 0, nullptr);
    SQLBindCol(stmt, 2, SQL_C_SLONG,  &instance.n_requests,       0, nullptr);
    SQLBindCol(stmt, 3, SQL_C_SLONG,  &instance.n_terminals,       0, nullptr);

    // Fetch first pickup row
    if (SQLFetch(stmt) != SQL_SUCCESS) {
        throw std::runtime_error("No instance information found");
    }

    /* ============================================================
    5. Cleanup ODBC handles
    ============================================================ */

    SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);

    return instance;
};


std::pair<double, double> compute_mrt(node_info pickup_node, node_info dropoff_node, double mrt_factor){
    double dist_pickup_dropoff = sqrt(pow(dropoff_node.x - pickup_node.x, 2) + pow(dropoff_node.y - pickup_node.y, 2));

    double mrt = dist_pickup_dropoff * mrt_factor;
    return {dist_pickup_dropoff, mrt};
    }


void update_time_windows(node_info& pickup_node, node_info& dropoff_node, double dist_pickup_dropoff ,double mrt){
    if (pickup_node.twL == 0){
        // std::cout << "Updating pickup time window // ";
        pickup_node.twU = std::ceil(dropoff_node.twU - dist_pickup_dropoff);
        pickup_node.twL = std::floor(dropoff_node.twL - mrt);
    }

    if (dropoff_node.twL == 0){
        // std::cout << "Updating dropoff time window // ";
        dropoff_node.twU = std::ceil(pickup_node.twU + mrt);
        dropoff_node.twL = std::floor(pickup_node.twL + dist_pickup_dropoff);
    }

}


std::pair<std::vector<node_info>, std::vector<node_info>> link_pickup_dropoff(const std::string& connStr, int instance_id, instance_info instance, double mrt_factor)
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
       2. Create ODBC connection
       ============================================================ */

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    // Allocate environment
    SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
    SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);

    // Allocate connection
    SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);

    // Connect using connection string
    SQLDriverConnectA(
        dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );


    /* ============================================================
       3. Fetch PICKUP node
       ============================================================ */

    for (int pickup_node_id = 0; pickup_node_id < instance.n_requests; ++pickup_node_id){

        node_info pickup{};

        SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
        SQLPrepareA(stmt, (SQLCHAR*)sql_pickup.c_str(), SQL_NTS);

        // Bind instance_id parameter
        SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &instance_id, 0, nullptr);

        SQLBindParameter(stmt, 2, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &pickup_node_id, 0, nullptr);

        // Execute query
        SQLExecute(stmt);

        // Bind result columns
        SQLBindCol(stmt, 1, SQL_C_SLONG,  &pickup.node_id, 0, nullptr);
        SQLBindCol(stmt, 2, SQL_C_DOUBLE, &pickup.x,       0, nullptr);
        SQLBindCol(stmt, 3, SQL_C_DOUBLE, &pickup.y,       0, nullptr);
        SQLBindCol(stmt, 4, SQL_C_DOUBLE, &pickup.twL,     0, nullptr);
        SQLBindCol(stmt, 5, SQL_C_DOUBLE, &pickup.twU,     0, nullptr);

        // Fetch first pickup row
        if (SQLFetch(stmt) != SQL_SUCCESS) {
            throw std::runtime_error("No pickup node found");
        }

        SQLFreeHandle(SQL_HANDLE_STMT, stmt);

    /* ============================================================
       4. Fetch DROPOFF node corresponding to pickup
       ============================================================ */

        node_info dropoff{};
        int dropoff_node_id = pickup_node_id + instance.n_requests;

        SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
        SQLPrepareA(stmt, (SQLCHAR*)sql_dropoff.c_str(), SQL_NTS);

        // Bind parameters
        SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &instance_id, 0, nullptr);

        SQLBindParameter(stmt, 2, SQL_PARAM_INPUT,
                        SQL_C_SLONG, SQL_INTEGER,
                        0, 0, &dropoff_node_id, 0, nullptr);

        // Execute query
        SQLExecute(stmt);

        // Bind result columns
        SQLBindCol(stmt, 1, SQL_C_SLONG,  &dropoff.node_id, 0, nullptr);
        SQLBindCol(stmt, 2, SQL_C_DOUBLE, &dropoff.x,       0, nullptr);
        SQLBindCol(stmt, 3, SQL_C_DOUBLE, &dropoff.y,       0, nullptr);
        SQLBindCol(stmt, 4, SQL_C_DOUBLE, &dropoff.twL,     0, nullptr);
        SQLBindCol(stmt, 5, SQL_C_DOUBLE, &dropoff.twU,     0, nullptr);

        if (SQLFetch(stmt) != SQL_SUCCESS) {
            throw std::runtime_error("No dropoff node found");
        }

        SQLFreeHandle(SQL_HANDLE_STMT, stmt);

        auto [dist_pickup_dropoff, mrt] = compute_mrt(pickup, dropoff, mrt_factor);

        std::cout << "Time windows prior to the update: [" << pickup.twL << "," << pickup.twU << "] // [" << dropoff.twL << "," << dropoff.twU << "] \n";
        update_time_windows(pickup, dropoff, dist_pickup_dropoff, mrt);
        std::cout << "Time windows after the update: [" << pickup.twL << "," << pickup.twU << "] // [" << dropoff.twL << "," << dropoff.twU << "] \n";


        pickup_nodes.push_back(pickup);
        dropoff_nodes.push_back(dropoff);

    };


    /* ============================================================
       5. Cleanup ODBC handles
       ============================================================ */

    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);




    /* ============================================================
       6. Return both nodes
       ============================================================ */

    return { pickup_nodes, dropoff_nodes };
}


void input_into_db(std::vector<node_info> pickup_nodes, std::vector<node_info> dropoff_nodes, std::string connStr, int instance_id){
    std::string sql_update_pickup = 
    "UPDATE dbo.nodes "
    "SET twL_constrained = ?, twU_constrained = ? "
    "WHERE instance_id = ? AND node_id = ?";

    std::string sql_update_dropoff = 
    "UPDATE dbo.nodes "
    "SET twL_constrained = ?, twU_constrained = ? "
    "WHERE instance_id = ? AND node_id = ?";

    /* ============================================================
    2. Create ODBC connection
    ============================================================ */

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    // Allocate environment
    SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
    SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);

    // Allocate connection
    SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);

    // Connect using connection string
    SQLDriverConnectA(
        dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );

    /* ============================================================
    5. Interact with the DB
    ============================================================ */

    // Update Pickup Windows

    for (node_info pickup_node : pickup_nodes){

    SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
    SQLPrepareA(stmt, (SQLCHAR*)sql_update_pickup.c_str(), SQL_NTS);

    // Bind parametres
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &pickup_node.twL, 0, nullptr);
    SQLBindParameter(stmt, 2, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &pickup_node.twU, 0, nullptr);
    SQLBindParameter(stmt, 3, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &instance_id, 0, nullptr);
    SQLBindParameter(stmt, 4, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &pickup_node.node_id, 0, nullptr);

    // Execute Query
    SQLExecute(stmt);

    SQLFreeHandle(SQL_HANDLE_STMT, stmt);

    }


    for (node_info dropoff_node : dropoff_nodes){

    // Update Dropoff Windows
    SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);
    SQLPrepareA(stmt, (SQLCHAR*)sql_update_dropoff.c_str(), SQL_NTS);

    // Bind parametres
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &dropoff_node.twL, 0, nullptr);
    SQLBindParameter(stmt, 2, SQL_PARAM_INPUT, SQL_C_DOUBLE, SQL_FLOAT, 0, 0, &dropoff_node.twU, 0, nullptr);
    SQLBindParameter(stmt, 3, SQL_PARAM_INPUT, SQL_C_SLONG, SQL_INTEGER, 0, 0, &dropoff_node.node_id, 0, nullptr);

    // Execute Query
    SQLExecute(stmt);

    SQLFreeHandle(SQL_HANDLE_STMT, stmt);}

    /* ============================================================
    5. Cleanup ODBC handles
    ============================================================ */

    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);
}


int main(){
    double mrt_factor = 5;

    std::string connStr = produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");
    add_column_to_table(connStr, "dbo.nodes", "twL_constrained", "FLOAT");
    add_column_to_table(connStr, "dbo.nodes", "twU_constrained", "FLOAT");

    std::vector<int> instance_ids = get_instances(connStr);

    // for (int instance_id : instance_ids){

    int instance_id = 1;

    std::cout << "Modifying DB for instance " << instance_id << "\n";

    instance_info instance = get_instance_node_info(connStr, instance_id);

    auto [pickup_nodes, dropoff_nodes] = link_pickup_dropoff(connStr, instance_id, instance, mrt_factor);

    input_into_db(pickup_nodes, dropoff_nodes, connStr, instance_id);

    std::cout << "Database modified for instance " << instance_id << "\n";

// }

    return 0;

}

