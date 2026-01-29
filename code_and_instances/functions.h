#include <string>
#include <iostream>

/* IDENTIFICATION OF THE NODES
------------------------------------------------------------------------------------------------------
0. Pickup nodes:			         0			    ->      n_requests - 1
1. Delivery nodes:		    n_requests			    ->  2 * n_requests - 1		  
2. Depot nodes:			2 * n_requests			    ->  2 * n_requests + n_vehicles - 1
3. Transfer nodes:		2 * n_requests + n_vehicles	->  2 * n_requests + n_vehicles + n_terminals - 1
------------------------------------------------------------------------------------------------------ */

/* IDENTIFICATION OF THE REFERENCES (for user i)
------------------------------------------------------------------------------------------------------
0. Pickup reference:				i
1. Delivery reference:				i +     n_requests
2. First transfer reference:		i + 2 * n_requests + 2 * n_vehicles
3. Second transfer reference:		i + 3 * n_requests + 2 * n_vehicles
------------------------------------------------------------------------------------------------------ */


// Structs related to the instance:
struct node { double x_coord; double y_coord; double lower_tw; double upper_tw; double service_dur; int type; };
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


// Preprocessing:
void build_problem(problem &p, std::string data_darp);
void read_data_darp(problem &p, const std::string& data_darp);
void compute_matrix_nontightened(problem &p);
void compute_matrix_pt(problem &p);
void compute_nearest_terminals(problem &p);
void compute_relative_mrt(problem &p);
void tighten_time_windows(problem &p);
void tighten_matrix(problem &p);
void compute_user_dissimilarities(problem &p);
void compute_node_dissimilarities(problem &p);

// Solution management:
void initialize_solution(problem &p, solution &s);
void update_solution(problem &p, solution &s1, solution &s2, bool create_schedule);

// Removal operators:
void random_removal(problem &p, solution &s, double removal_percentage, bool clustered_removal);
void worst_removal(problem &p, solution &s, double removal_percentage, bool clustered_removal, bool idarp_adapted);
void related_removal(problem &p, solution &s, double removal_percentage, bool clustered_removal, bool idarp_adapted);
void route_removal(problem &p, solution &s, double removal_percentage, bool clustered_removal);

// Insertion operators:
void random_order_insertion(problem &p, solution &s, double max_objective_value = std::numeric_limits<double>::infinity());
void greedy_insertion(problem &p, solution &s, double max_objective_value = std::numeric_limits<double>::infinity());
void two_regret_insertion(problem &p, solution &s, double max_objective_value = std::numeric_limits<double>::infinity());

// Local search operators:
void relocate(problem &p, solution &s);
void exchange_natural_sequences(problem &p, solution &s);

// Supporting functions:
void remove_request(problem &p, solution &s, int request);
int clear_terminal(problem &p, solution &s, int request);

// Supporting functions for insertions:
void clear_insertion_data(problem &p, solution &s);
void construct_request_banks(problem &p, solution &s);
void compute_best_insertion_request(problem &p, solution &s, int request, double max_objective_value, bool track_second_best = false);
void update_best_insertion_request(struct problem &p, struct solution &s, int request, int update_vehicle1, int update_vehicle2, double max_objective_value, bool track_second_best = false);
void compute_best_insertion_without_pt(problem &p, solution &s, int request, int vehicle, bool track_second_best);
void compute_best_insertion_pt(problem &p, solution &s, int request, int vehicle1, int vehicle2, int terminal1, int termninal2, bool track_second_best);
void insert_request(problem &p, solution &s, int request, int vehicle);
void insert_request_with_transfer(problem &p, solution &s, int request, int vehicle1, int vehicle2);

// Feasibility checks:
bool check_load(problem &p, solution &s, int ref_first_node, int ref_last_node);
bool check_completeness(problem &p, solution &s);
bool check_service_start_Tang(problem &p, solution &s, int vehicle, bool update = false);
bool check_service_start_BellmanFord(problem &p, solution &s, bool update = false);
double compute_required_time_pt(problem &p, int node_origin_terminal, int node_destination_terminal, double pickup_time_second_transfer);

// Update solution attributes:
void update_load(problem &p, solution &s, int ref_first_node, int ref_last_node);
void update_earliest_time(problem &p, solution &s, int begin_node);
void update_latest_time(problem &p, solution &s, int end_node);

// Additional check to prove the feasibility:
bool check_feasibility(problem &p, solution &s);

// Write output file:
void write_output_file(problem &p, solution &s, int iterations, double computation_time, bool extensive_output, std::string file);
void write_solution_to_file(const std::string& file, const solution& s, int n_nodes);