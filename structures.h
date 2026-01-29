#pragma once

// Necessary Libraries
#include <string>
#include <windows.h>
#include <sqlext.h>
#include <vector>



struct sql_params {
    SQLHENV  env  = SQL_NULL_HENV;
    SQLHDBC  dbc  = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;
};

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

struct instance_info{
    int instance_id;
    int n_vehicles;
    int n_requests;
    int n_terminals;
};

struct Train {
            int id;
            std::string stations;
            double departure;
            double arrival;
        };

struct request_info{
    int request_id;
    node_info pickup_node;
    node_info dropoff_node;
};

struct db_request_info{
    int request_id;
    int instance_id;
    int pickup_node_id;
    int dropoff_node_id;
    double t_pickup_dropoff;
    double mrt;
    int passengers;
    double twL_pickup;
    double twU_pickup;
    double twL_dropoff;
    double twU_dropoff;
};

// Structs related to the instance:
struct node { double x_coord; double y_coord; double lower_tw; double upper_tw; double service_dur; int type; int actual_node; }; //double service_start_time;};
struct request { int users_mobile; int users_wheelchair; double max_ride_time; bool eligible_pt; bool long_dist; int nearest_t1[5]; int nearest_t2[5]; };
struct vehicle { int cap_mobile; int cap_wheelchair; double max_route_dur; };

// Struct defining the problem:
struct problem { int n_requests; int n_vehicles; int n_terminals; int n_nodes; int n_references;
				 double *matrix; double *matrix_nontightened; double *user_dissim; double *node_dissim; double *matrix_pt;
				 node *nodes; request *requests; vehicle *vehicles; };

// Struct defining the solution:
struct solution { double total_distance; int *predecessor; int *successor;
				  int *map_to_node; int *map_to_request; int *map_to_vehicle;
				  int *load_mobile; int *load_wheelchair;
				  double *service_start; double *earliest_time; double *latest_time;
				  double *best_insertion_cost; double *second_best_insertion_cost;
				  int *best_predecessor; int *best_successor;
				  int *best_vehicle1; int *best_vehicle2; int *second_best_vehicle1; int *second_best_vehicle2;
				  int *best_terminal1; int *best_terminal2;
				  std::vector<int> request_bank_prior; std::vector<int> request_bank_non_prior;
};


struct removable_node_info {
    int reference;
    int node_id_problem;
    int vehicle_number;
    int request_id;
    std::string node_type;
};
