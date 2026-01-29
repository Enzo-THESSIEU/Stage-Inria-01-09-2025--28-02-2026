// Inclusions
#include "Pipelines.h"
#include "SQL_Queries.h"
#include "information_modifying_functions.h"
#include "db_connection_functions.h"
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
#include <unordered_map>
#include <cstdlib>
#include <algorithm>

///g++ -std=c++17 -O2 -Wall -Wextra \
Pipelines.cpp \
SQL_Queries.cpp \
information_modifying_functions.cpp \
db_connection_functions.cpp \
-o Pipelines.exe -lodbc32 -lodbccp32

/// ./Pipelines.exe



void import_results_pipeline(){

    odbc_functions odbc;
    table_creators creators;
    input_to_db importer;

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

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    creators.create_results_table_if_missing(sql_param);

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    int instance_id = 1;

    for (const std::string& file : data_darp){
        std::string path = "C:\\Users\\enzot\\Documents\\Cesure\\1ere cesure inria Lille\\Codes\\Stage-Inria-01-09-2025--28-02-2026\\code_and_instances\\" + file ;
        importer.import_results_summary_tsv(path, instance_id, sql_param);
        instance_id++ ;
    };
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
    sql_param = odbc.cleanup_sql_query(sql_param);

};



void import_most_correlated_requests_pipiline(){
    odbc_functions odbc;
    table_creators creators;
    input_to_db importer;
    info_getters getter;
    functions function;

    double gamma_WT = 1; 
    double gamma_TW = 1;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    creators.create_request_table_if_missing(sql_param);

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    std::vector<int> instances = getter.get_instances(sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    for (int instance : instances){
        instance_info info_of_instance = getter.get_instance_info(instance, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        int inst_id = info_of_instance.instance_id;
        int n_requests = info_of_instance.n_requests;
        std::vector<request_info> all_requests;
        for (int request_id = 0 ; request_id < n_requests ; request_id++){
            request_info info_of_request;
            info_of_request.request_id = request_id;
            info_of_request.pickup_node = getter.get_node_info(request_id, inst_id, sql_param);
            sql_param = odbc.free_up_statement_for_new_query_type(sql_param);       
            info_of_request.dropoff_node = getter.get_node_info(request_id + n_requests, inst_id, sql_param);
            sql_param = odbc.free_up_statement_for_new_query_type(sql_param);   
            all_requests.push_back(info_of_request);
        };
        for (request_info request : all_requests){
            std::vector<int> correlated_set;
            correlated_set = function.correlation_set_request(request, gamma_WT, gamma_TW, all_requests);
            importer.insert_correlated_set(sql_param, correlated_set, inst_id, request.request_id);
            sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        };
    };
    sql_param = odbc.cleanup_sql_query(sql_param);
};



void update_TW_pipeline(){
    odbc_functions odbc;
    table_creators creators;
    input_to_db importer;
    info_getters getter;
    functions function;

    double mrt_factor = 1.5;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    std::vector<int> instances = getter.get_instances(sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
    for (int instance : instances){
        // int instance = 1 ;
        instance_info info_of_instance = getter.get_instance_info(instance, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        int n_nodes = info_of_instance.n_requests +  info_of_instance.n_terminals + info_of_instance.n_vehicles;
        std::vector<node_info> pickup_nodes = {};
        std::vector<node_info> dropoff_nodes = {};
        for (int i = 0; i < info_of_instance.n_requests; i++){
            node_info pickup_node = getter.get_node_info(instance, i, sql_param);
            sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
            node_info dropoff_node = getter.get_node_info(instance, i + info_of_instance.n_requests, sql_param);
            sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
            auto [dist_pickup_dropoff, mrt] = function.compute_mrt(pickup_node, dropoff_node, mrt_factor);
            function.update_time_windows(pickup_node, dropoff_node, dist_pickup_dropoff ,mrt);

            // std::cout << "Time windows for request " << pickup_node.node_id << "\n";
            // std::cout << "tij = " << dist_pickup_dropoff << ", mrt = " << mrt << "\n";
            // std::cout << "Pickup node [" << pickup_node.twL << "," << pickup_node.twU << "] // Dropoff node [" << dropoff_node.twL << "," << dropoff_node.twU << "] \n";
            // std::cout << "Pickup node [" << pickup_node.twL_constrained << "," << pickup_node.twU_constrained << "] // Dropoff node [" << dropoff_node.twL_constrained << "," << dropoff_node.twU_constrained << "] \n";
            // std::cout << "------------------------------------------------------------------------------------------- \n";

            pickup_nodes.push_back(pickup_node);
            dropoff_nodes.push_back(dropoff_node);
        };
        importer.input_constrained_TW_into_db(pickup_nodes, dropoff_nodes, instance, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
    };
    sql_param = odbc.cleanup_sql_query(sql_param);
};



void build_TW_score_for_subset(){
    odbc_functions odbc;
    table_creators creators;
    table_modifiers modifier;
    input_to_db importer;
    info_getters getter;
    functions function;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    creators.create_request_table_if_missing(sql_param);

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    int inst_id = 1;

    instance_info info_of_instance = getter.get_instance_info(inst_id, sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
    int n_requests = info_of_instance.n_requests;
    std::vector<request_info> all_requests;
    std::vector<std::pair<int, double>> tw_scores;
    for (int request_id = 0 ; request_id < n_requests ; request_id++){
        request_info info_of_request;
        info_of_request.request_id = request_id;
        info_of_request.pickup_node = getter.get_node_info(inst_id, request_id, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);       
        info_of_request.dropoff_node = getter.get_node_info(inst_id, request_id + n_requests, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);   
        all_requests.push_back(info_of_request);
        double tw_score = function.time_window_score(info_of_request);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        importer.input_to_db_tw_scores(sql_param, tw_score, inst_id, request_id);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param); 
    };

    odbc.cleanup_sql_query(sql_param);



};



void build_request_db_pipeline(){
    odbc_functions odbc;
    table_creators creators;
    table_modifiers modifier;
    input_to_db importer;
    info_getters getter;
    functions function;

    double mrt_factor = 1.5;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    creators.delete_table(sql_param, "dbo.requests");

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    creators.create_request_table_if_missing(sql_param);

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    int inst_id = 1;

    instance_info info_of_instance = getter.get_instance_info(inst_id, sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
    int n_requests = info_of_instance.n_requests;
    request_info req_info;
    for (int request_id = 0 ; request_id < n_requests ; request_id++){
        request_info info_of_request;
        info_of_request.request_id = request_id;
        info_of_request.pickup_node = getter.get_node_info(inst_id, request_id, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);       
        info_of_request.dropoff_node = getter.get_node_info(inst_id, request_id + n_requests, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);  
        
        db_request_info db_info_of_req;
        db_info_of_req.request_id = request_id;
        db_info_of_req.instance_id = inst_id;
        db_info_of_req.pickup_node_id = info_of_request.pickup_node.node_id;
        db_info_of_req.dropoff_node_id = info_of_request.dropoff_node.node_id;
        auto[t_p_d, mrt] = function.compute_mrt(info_of_request.pickup_node, info_of_request.dropoff_node, mrt_factor);
        db_info_of_req.t_pickup_dropoff = t_p_d;
        db_info_of_req.mrt = mrt;
        db_info_of_req.passengers = info_of_request.pickup_node.mobile + info_of_request.pickup_node.wheelchair;
        db_info_of_req.twL_pickup = info_of_request.pickup_node.twL_constrained;
        db_info_of_req.twU_pickup = info_of_request.pickup_node.twU_constrained;
        db_info_of_req.twL_dropoff = info_of_request.dropoff_node.twL_constrained;
        db_info_of_req.twU_dropoff = info_of_request.dropoff_node.twU_constrained;

        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        importer.input_to_db_request_info(sql_param, db_info_of_req);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param); 
    };
    odbc.cleanup_sql_query(sql_param);
};



std::string print_n_tw_score_requests(int instance_id, int subset_number, int size_of_subset, int size_of_set, int offset){
    odbc_functions odbc;
    table_creators creators;
    table_modifiers modifier;
    input_to_db importer;
    info_getters getter;
    functions function;

    double mrt_factor = 1.5;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

     

    std::vector<db_request_info> requests = getter.get_n_requests(size_of_subset, offset, instance_id, sql_param);

    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    std::vector<node_info> pickups;  // node + passengers
    std::vector<node_info> dropoffs; // node + passengers

    for (const db_request_info& request : requests) {
        node_info pickup_node = getter.get_node_info(request.instance_id, request.pickup_node_id, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        node_info dropoff_node = getter.get_node_info(request.instance_id, request.dropoff_node_id, sql_param);
        sql_param = odbc.free_up_statement_for_new_query_type(sql_param);
        pickups.push_back(pickup_node);
        dropoffs.push_back(dropoff_node);
    };

    std::vector<node_info> depots    = getter.get_nodes_by_type(instance_id, "Depot Node",    21, sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    std::vector<node_info> terminals = getter.get_nodes_by_type(instance_id, "Transfer Node", 81, sql_param);
    sql_param = odbc.free_up_statement_for_new_query_type(sql_param);

    // for (const db_request_info& request : requests) {
    //     std::cout
    //         << "Request id: " << request.request_id
    //         << " | instance: " << request.instance_id
    //         << " | pickup_node: " << request.pickup_node_id
    //         << " | dropoff_node: " << request.dropoff_node_id
    //         << " | t_pd: " << request.t_pickup_dropoff
    //         << " | mrt: " << request.mrt
    //         << " | passengers: " << request.passengers
    //         << " | tw_pickup: [" << request.twL_pickup << ", " << request.twU_pickup << "]"
    //         << " | tw_dropoff: [" << request.twL_dropoff << ", " << request.twU_dropoff << "]"
    //         << '\n';
    // };

    std::vector<node_info> nodes;
    nodes.reserve(pickups.size() + dropoffs.size() + depots.size() + terminals.size());

    nodes.insert(nodes.end(), pickups.begin(), pickups.end());
    nodes.insert(nodes.end(), dropoffs.begin(), dropoffs.end());
    nodes.insert(nodes.end(), depots.begin(),  depots.end());
    nodes.insert(nodes.end(), terminals.begin(), terminals.end());

    std::string file_name = "datafile_subset_" + std::to_string(subset_number) + ".txt";

    std::ofstream out(file_name, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        odbc.cleanup_sql_query(sql_param);
        throw std::runtime_error("Could not create datafile_subset_" + std::to_string(subset_number) + ".txt");
    }

    out << depots.size() << " "
        << size_of_subset << " "
        << terminals.size() << " "
        << "\n";

    auto write_node = [&](const node_info& node, int k){
        out << k << " "
            << node.x << " " << node.y << " "
            << node.twL_constrained << " " << node.twU_constrained << " "
            << node.mobile << " "
            << node.wheelchair << " "
            << node.service_dur << " "
            << node.node_id << " "
            << "\n";
    };

    auto write_non_request_node = [&](const node_info& node, int k){
        out << k << " "
            << node.x << " " << node.y << " "
            << node.twL << " " << node.twU << " "
            << node.mobile << " "
            << node.wheelchair << " "
            << node.service_dur << " "
            << node.node_id << " "
            << "\n";
    };

    int max_k = nodes.size();
    for (int k = 0 ; k < max_k ; k++){
        if (k < 2 * size_of_subset) {
            // std::cout << "Writing Request Nodes for k = " << k << " < " << 2 * size_of_subset << "\n";
            write_node(nodes[k], k);}
        else {
            // std::cout << "Writing non Request Nodes for k = " << k << " >= " << 2 * size_of_subset << "\n";
            write_non_request_node(nodes[k], k);}
    }

    return file_name;

};





instance_info get_instance_info_pipeline(int instance_id){
    odbc_functions odbc;
    table_creators creators;
    input_to_db importer;
    info_getters getter;
    functions function;

    double mrt_factor = 1.5;

    std::string connStr = odbc.produce_conn_str(R"(localhost\SQLEXPRESS)", "darp");

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    sql_params sql_param = {env, dbc, stmt};

    sql_param = odbc.connect_to_db(connStr, sql_param);

    instance_info info_of_instance = getter.get_instance_info(instance_id, sql_param);

    odbc.free_up_statement_for_new_query_type(sql_param);
    odbc.cleanup_sql_query(sql_param);

    return info_of_instance;
};

// ------------------ MAIN ------------------
// int main(){
//     // update_TW_pipeline();
//     build_request_db_pipeline();
//     build_TW_score_for_subset();
//     print_n_tw_score_requests();
//     return 0;
// };






























































