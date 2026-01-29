#include <windows.h>
#include <sqlext.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <cstdlib>


// ---------------- Structure Definition -------------------
struct node_info{
    int node_id;
    double x;
    double y;
    double twL;
    double twU;
    int mobile;
    int wheelchair;
    double service_dur;
    std::string node_type;
    double twL_constrained;
    double twU_constrained;
};

struct sql_params{
    SQLHENV env;
    SQLHDBC dbc;
    SQLHSTMT stmt; 
};

struct instance_info{
    int id;
    int n_vehicles;
    int n_requests;
    int n_terminals;
};

struct additional_node_info{
    double x;
    double y;

    double tw_mid;
    double tw_width;

    double dist_to_nearest_depot;
    double dist_to_nearest_transfer;
};

struct additional_request_info{
    double px;
    double py;
    double dx;
    double dy;

    double twp_mid;
    double twd_mid;

    double twp_width;
    double twd_width;

    double dist_p_d;
    double dist_to_nearest_depot;
    double dist_to_nearest_transfer;
};

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
    SQLDriverConnectA(
        dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );

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


// ----------------- SQL Query Functions ---------------------
node_info get_node_info(int instance_id, int node_id, sql_params sql_param){
    node_info node{};

    SQLHENV env = sql_param.env;
    SQLHDBC dbc = sql_param.dbc;
    SQLHSTMT stmt = sql_param.stmt;

    std::string sql = "SELECT node_id, x, y, twL, twU, mobile, wheelchair, service_dur, node_type, twL_constrained, twU_constrained "
                        "FROM dbo.nodes "
                        "WHERE instance_id = ? AND node_id =?;";

    SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Bind Parametres

    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                    SQL_C_SLONG, SQL_INTEGER,
                    0, 0, &instance_id, 0, nullptr);

    SQLBindParameter(stmt, 2, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &node_id, 0, nullptr);

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

        // Execute query
        SQLExecute(stmt);

        // Fetch

        if (SQLFetch(stmt) != SQL_SUCCESS) {
            throw std::runtime_error("Node {node_id} not found");
        };

    node.node_type = node_type_buf;

    return node;


};


instance_info get_instance_info(int instance_id, sql_params sql_param){
    instance_info instance{};

    SQLHENV env = sql_param.env;
    SQLHDBC dbc = sql_param.dbc;
    SQLHSTMT stmt = sql_param.stmt;

    std::string sql = 
    "SELECT id, n_vehicles, n_requests, n_terminals "
    "FROM dbo.instances "
    "WHERE id = ?";

    SQLPrepareA(stmt, (SQLCHAR*)sql.c_str(), SQL_NTS);

    // Bind instance_id parameter
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT,
                     SQL_C_SLONG, SQL_INTEGER,
                     0, 0, &instance_id, 0, nullptr);

    // Bind result columns
    SQLBindCol(stmt, 1, SQL_C_SLONG,  &instance.id, 0, nullptr);
    SQLBindCol(stmt, 2, SQL_C_SLONG,  &instance.n_vehicles, 0, nullptr);
    SQLBindCol(stmt, 3, SQL_C_SLONG,  &instance.n_requests,       0, nullptr);
    SQLBindCol(stmt, 4, SQL_C_SLONG,  &instance.n_terminals,       0, nullptr);

    // Execute query
    SQLExecute(stmt);

    // Fetch first pickup row
    if (SQLFetch(stmt) != SQL_SUCCESS) {
        throw std::runtime_error("No instance information found");
    };

    return instance;

};



// ----------------- Info Manipulation Functions -------------
double dij(node_info node_i, node_info node_j){
    return sqrt(pow(node_i.x - node_j.x, 2) + pow(node_i.y - node_j.y, 2));
}


double correlation_function(node_info node_i, node_info node_j, double gamma_WT, double gamma_TW){

    double dist_i_j = dij(node_i, node_j);

    double max_1 = std::max(node_j.twL_constrained - node_i.service_dur - dist_i_j - node_i.twU_constrained,  0.0);
    double max_2 = std::max(node_i.twL_constrained - node_i.service_dur - dist_i_j - node_j.twU_constrained,  0.0);

    return dist_i_j + gamma_WT * max_1 + gamma_TW * max_2;
};



bool same_request(node_info node_i, node_info node_j, instance_info instance){
    const bool i_isPickup   = (node_i.node_type == "pickup");
    const bool j_isPickup   = (node_j.node_type == "pickup");
    const bool i_isDropoff  = (node_i.node_type == "dropoff");
    const bool j_isDropoff  = (node_j.node_type == "dropoff");

    if ((i_isDropoff && j_isPickup) || (j_isDropoff && i_isPickup)){
        return std::abs(node_i.node_id - node_j.node_id) == instance.n_requests;
    };

    return false;
};


std::unordered_map<int,int> build_pair_pi_di(instance_info instance){
    std::unordered_map<int,int> pair_pi_di = {};

    for (int i = 0; i < instance.n_requests; i++){
         pair_pi_di[i] = i + instance.n_requests;
    };

    return pair_pi_di;
}


/*
Implements the same *arc eliminations* you had in Python, but as a boolean test:
return true  => arc i -> j is allowed
return false => arc i -> j is forbidden (would have been deleted from tij)

Inputs:
- node_i, node_j: the two nodes
- zeroDepotId / endDepotId: node_id of start depot and end depot
- isPickup/isDropoff/isDepot/isTransfer determined via node_type strings
- pair_pi_di: mapping pickup_id -> dropoff_id (regular requests)
- pair_pi_di_M: mapping pickup_id -> dropoff_id (MoPS, optional)
- Lbar: max ride time per pickup_id (keyed by pickup node_id)
- use_MoPS: enable MoPS-related rules
- do_timewindow_test / do_ride_time_test: those blocks were commented in Python; keep off by default unless you want them on
*/


bool viable_trip(
    const node_info& node_i,
    const node_info& node_j
){
    const int bi = node_i.node_id;
    const int bj = node_j.node_id;

    const bool i_isDepot    = (node_i.node_type == "depot");
    const bool j_isDepot    = (node_j.node_type == "depot");
    const bool i_isPickup   = (node_i.node_type == "pickup");
    const bool j_isPickup   = (node_j.node_type == "pickup");
    const bool i_isDropoff  = (node_i.node_type == "dropoff");
    const bool j_isDropoff  = (node_j.node_type == "dropoff");
    const bool i_isTransfer = (node_i.node_type == "transfer");
    const bool j_isTransfer = (node_j.node_type == "transfer");

    const double tij = dij(node_i, node_j);
    const double ei  = node_i.twL_constrained;   
    const double di  = node_i.service_dur;
    const double lj  = node_j.twU_constrained;   

    // ------------------------------------------------------------
    // 1) Depot restrictions  (your Python "Depot restrictions" block)
    //    - Remove arcs to start depot
    //    - Remove arcs from end depot
    //    - Remove start_depot -> any dropoff
    //    - Remove any pickup -> end_depot
    // ------------------------------------------------------------

    if (i_isDepot && j_isDropoff) return false; // (zeroDepot, D) deleted
    if (j_isDepot  && i_isPickup)  return false; // (P, endDepot) deleted

    // (Optional) self-loops were commented out in Python
    if (bi == bj) return false;

    // ------------------------------------------------------------
    // 3) Time-window infeasibility (COMMENTED OUT in your Python)
    //    Python condition: ei[i] + di[i] + tij(i,j) >= li[j] => delete(i,j)
    // ------------------------------------------------------------
  
    if (ei + di + tij >= lj) return false;
    
    // ------------------------------------------------------------
    // 4) Ride time infeasibility (COMMENTED OUT in your Python)
    //    Python: for i in P, j in N:
    //      if tij(i,j) + di[j] + tij(j, drop(i)) >= Lbar[i] => delete (i,j) and (j,drop(i))
    //    Here we can only decide viability of (i->j), so we just test the condition and return false.
    //    We need drop(i) => pair_pi_di[i]
    // ------------------------------------------------------------
    // if (i_isPickup){
    //     int drop_i = pair_pi_di.find(bi);
    //     if (tij + node_j.service_dur + dij(bj, drop_i) > mrt_factor * dij(bi, drop_i)){

    //     }
    // }

    return true; // arc kept
}


std::vector<std::pair<int, double>> correlation_set(node_info node_i, double gamma_WT, double gamma_TW, std::vector<node_info> nodes){
    std::vector<std::pair<int, double>> corr;

    for (const auto& node_j : nodes) {
        if (viable_trip(node_i, node_j)) { // or viable_trip(node_i, node_j)
            double score = correlation_function(node_i, node_j, gamma_WT, gamma_TW);
            corr.emplace_back(node_j.node_id, score);
        }
    }

    std::sort(corr.begin(), corr.end(), [](const auto& a, const auto& b){ return a.second < b.second; });

    return corr;
}





int main(){
    std::string connStr = produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = connect_to_db(connStr, sql_param);

    int instance_id = 1;
    int node_id = 1;
    auto node = get_node_info(instance_id, node_id, sql_param);
    sql_param = free_up_statement_for_new_query_type(sql_param);
    sql_param = cleanup_sql_query(sql_param);

};