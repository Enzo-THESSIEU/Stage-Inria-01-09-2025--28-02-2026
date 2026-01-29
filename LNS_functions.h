#include <string>
#include <iostream>
#include <vector>
#include <limits>
#include <utility>


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


// Structures
#include "structures.h"


// Preprocessing:
// void build_problem(problem &p, std::string data_darp);
void build_first_problem(struct problem &p, std::string data_darp);
void build_sequential_problem(struct problem &current_p, std::string data_darp, problem previous_p, solution &previous_s, std::vector<std::pair<int, int>> phantom_pickup_nodes, std::vector<std::pair<int, int>> phantom_dropoff_nodes, bool is_final_subset);

void read_data_darp(problem &p, const std::string& data_darp);
void read_sequential_data_darp(problem &current_p, const std::string& data_darp, solution previous_s, problem previous_p);

void compute_matrix_nontightened(problem &p);
void compute_matrix_pt(problem &p);
void compute_nearest_terminals(problem &p);
void compute_relative_mrt(problem &p);
void tighten_time_windows(problem &p);
void tighten_matrix(problem &p);
void compute_user_dissimilarities(problem &p);
void compute_node_dissimilarities(problem &p);

void compute_phantom_end_depot_node(struct problem &p, int phantom_node_id);
void compute_phantom_node_as_pickup_node(struct problem &p, int phantom_node_id, int entrynode);
void compute_phantom_node_dropoff_node(struct problem &p, int phantom_node_id, int exitnode);

// Subproblem Management
std::vector<int> get_intersecting_requests(std::string previous_file, int overlap_size);
std::vector<std::pair<int, int>> get_vehicle_routes(struct problem &p, struct solution &s, int vehicle_number);
std::pair<std::vector<std::pair<int,int>>, std::vector<std::pair<int,int>>> break_up_vehicle_route(solution s, problem &p, std::vector<std::pair<int,int>> vehicle_route);
static inline void unlink_ref(solution &s, int ref);
static inline removable_node_info *find_node_by_req(std::vector<removable_node_info> &v, int request_id);
std::pair<std::vector<std::pair<int, int>> , std::vector<std::pair<int, int>>> remove_nodes_route(solution &s, problem &p, std::string previous_file, std::string current_file, int overlap_size);
std::vector<std::pair<int, double>> get_final_vehicle_positions(solution s, problem &p);

// Subproblem reuniting
std::vector<std::vector<std::pair<int, double>>> reconstruct_full_routes(const std::vector<std::pair<problem, solution>> &all_subset_info);


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
void write_solution_to_file(const solution& s, int n_nodes, int subproblem_number, int previous);

std::vector<std::pair<int,int>> get_vehicle_routes(problem &p, solution &s, int vehicle_num);
void write_subproblem_solution_output_file(problem &p, solution &s,
                                 std::string file, 
								 int subproblem_number);

void write_subproblem_output_file(problem &p,
                                 std::string file, 
								 int subproblem_number);

void write_full_solution_results_file(const std::vector<std::pair<problem, solution>> &all_subset_info);

// Runner
solution run_for_subproblem(int argc, char *argv[], std::string file, problem p, int subproblem_number); 

void init_output_directory();