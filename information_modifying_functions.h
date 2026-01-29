#pragma once

// Necessary Libraries
#include "db_connection_functions.h"
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


// // Structures
// struct sql_params {
//     SQLHENV  env  = SQL_NULL_HENV;
//     SQLHDBC  dbc  = SQL_NULL_HDBC;
//     SQLHSTMT stmt = SQL_NULL_HSTMT;
// };

// struct node_info{
//     int node_id;
//     double x;
//     double y;
//     double twL;
//     double twU;
//     int mobile;
//     int wheelchair;
//     double service_dur;
//     std::string node_type;
//     double twL_constrained;
//     double twU_constrained;
// };

// struct instance_info{
//     int instance_id;
//     int n_vehicles;
//     int n_requests;
//     int n_terminals;
// };

// struct Train {
//             int id;
//             std::string stations;
//             double departure;
//             double arrival;
//         };

// struct request_info{
//     int request_id;
//     node_info pickup_node;
//     node_info dropoff_node;
// };



class functions{
public:
    static double dij(node_info node_i, node_info node_j);

    std::pair<double, double> compute_mrt(node_info pickup_node, node_info dropoff_node, double mrt_factor);
    void update_time_windows(node_info& pickup_node, node_info& dropoff_node, double dist_pickup_dropoff ,double mrt);
    // std::string minutes_to_datetime(double minutes, const std::string& date = "2026-01-01");

    bool is_same_request(node_info node_i, node_info node_j, instance_info instance);
    // std::unordered_map<int,int> build_pair_pi_di(instance_info instance);
    bool viable_trip(const node_info& node_i,const node_info& node_j);

    double correlation_function(node_info node_i, node_info node_j, double gamma_WT, double gamma_TW);
    std::vector<std::pair<int, double>> correlation_set_node(node_info node_i, double gamma_WT, double gamma_TW, std::vector<node_info> nodes);
    std::vector<int> correlation_set_request(request_info request, double gamma_WT, double gamma_TW, std::vector<request_info> all_requests);
    double time_window_score(request_info request);
};