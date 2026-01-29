// Inclusions
#include "information_modifying_functions.h"

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

// Functions
double functions::dij(node_info node_i, node_info node_j){
    return sqrt(pow(node_i.x - node_j.x, 2) + pow(node_i.y - node_j.y, 2));
}



bool functions::is_same_request(node_info node_i, node_info node_j, instance_info instance){
    const bool i_isPickup   = (node_i.node_type == "pickup");
    const bool j_isPickup   = (node_j.node_type == "pickup");
    const bool i_isDropoff  = (node_i.node_type == "dropoff");
    const bool j_isDropoff  = (node_j.node_type == "dropoff");

    if ((i_isDropoff && j_isPickup) || (j_isDropoff && i_isPickup)){
        return std::abs(node_i.node_id - node_j.node_id) == instance.n_requests;
    };

    return false;
};



// std::unordered_map<int,int> functions::build_pair_pi_di(instance_info instance){
//     std::unordered_map<int,int> pair_pi_di = {};
//     for (int i = 0; i < instance.n_requests; i++){
//          pair_pi_di[i] = i + instance.n_requests;
//     };
//     return pair_pi_di;
// };



bool functions::viable_trip(const node_info& node_i, const node_info& node_j){
    const int bi = node_i.node_id;
    const int bj = node_j.node_id;

    const bool i_isDepot    = (node_i.node_type == "depot");
    const bool j_isDepot    = (node_j.node_type == "depot");
    const bool i_isPickup   = (node_i.node_type == "pickup");
    const bool j_isDropoff  = (node_j.node_type == "dropoff");

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



double functions::correlation_function(node_info node_i, node_info node_j, double gamma_WT, double gamma_TW){

    double dist_i_j = dij(node_i, node_j);

    double max_1 = std::max(node_j.twL_constrained - node_i.service_dur - dist_i_j - node_i.twU_constrained,  0.0);
    double max_2 = std::max(node_i.twL_constrained - node_i.service_dur - dist_i_j - node_j.twU_constrained,  0.0);

    return dist_i_j + gamma_WT * max_1 + gamma_TW * max_2;
};

std::vector<std::pair<int, double>> functions::correlation_set_node(node_info node_i, double gamma_WT, double gamma_TW, std::vector<node_info> nodes){
    std::vector<std::pair<int, double>> corr;

    for (const auto& node_j : nodes) {
        if (viable_trip(node_i, node_j)) { // or viable_trip(node_i, node_j)
            double score = correlation_function(node_i, node_j, gamma_WT, gamma_TW);
            corr.emplace_back(node_j.node_id, score);
        }
    }

    

    return corr;
}

std::vector<int> functions::correlation_set_request(request_info request, double gamma_WT, double gamma_TW, std::vector<request_info> all_requests){
    node_info pickup_node = request.pickup_node;
    node_info dropoff_node = request.dropoff_node;

    std::vector<std::pair<int, double>> corr_request;

    for (request_info other_request : all_requests){
        if (other_request.request_id == request.request_id){
            continue;
        }
        node_info other_request_pickup = other_request.pickup_node;
        node_info other_request_dropoff = other_request.dropoff_node;
        double score = correlation_function(pickup_node, other_request_pickup, gamma_WT, gamma_TW);
        score += correlation_function(pickup_node, other_request_dropoff, gamma_WT, gamma_TW);
        score += correlation_function(dropoff_node, other_request_pickup, gamma_WT, gamma_TW);
        score += correlation_function(dropoff_node, other_request_dropoff, gamma_WT, gamma_TW);
        corr_request.emplace_back(other_request.request_id, score);
    };

    std::sort(corr_request.begin(), corr_request.end(), [](const auto& a, const auto& b){ return a.second < b.second; });

        // Return only node_ids, ordered
    std::vector<int> result;
    result.reserve(corr_request.size());
    for (const auto& [id, _] : corr_request)
        result.push_back(id);

    return result;

}

double functions::time_window_score(request_info request){
    return (request.pickup_node.twU_constrained + request.pickup_node.twL_constrained + request.dropoff_node.twU_constrained + request.dropoff_node.twL_constrained) / 4;
};



std::pair<double, double> functions::compute_mrt(node_info pickup_node, node_info dropoff_node, double mrt_factor){
    double dist_pickup_dropoff = sqrt(pow(dropoff_node.x - pickup_node.x, 2) + pow(dropoff_node.y - pickup_node.y, 2));

    double mrt = dist_pickup_dropoff * mrt_factor;
    return {dist_pickup_dropoff, mrt};
    }

void functions::update_time_windows(node_info& pickup_node, node_info& dropoff_node, double dist_pickup_dropoff ,double mrt){
    if (pickup_node.twL == 0){
        // std::cout << "Updating pickup time window // ";
        pickup_node.twU_constrained = std::ceil(dropoff_node.twU - dist_pickup_dropoff);
        pickup_node.twL_constrained = std::floor(dropoff_node.twL - mrt);
        dropoff_node.twU_constrained = dropoff_node.twU;
        dropoff_node.twL_constrained = dropoff_node.twL;
    }

    if (dropoff_node.twL == 0){
        // std::cout << "Updating dropoff time window // ";
        pickup_node.twU_constrained = pickup_node.twU;
        pickup_node.twL_constrained = pickup_node.twL;
        dropoff_node.twU_constrained = std::ceil(pickup_node.twU + mrt);
        dropoff_node.twL_constrained = std::floor(pickup_node.twL + dist_pickup_dropoff);
    }

}



// std::string functions::minutes_to_datetime(double minutes, const std::string& date = "2026-01-01") {
//     int total = (int)minutes;
//     int hh = total / 60;
//     int mm = total % 60;
//     char buf[32];
//     std::snprintf(buf, sizeof(buf), "%s %02d:%02d:00", date.c_str(), hh, mm);
//     return std::string(buf);
// }




