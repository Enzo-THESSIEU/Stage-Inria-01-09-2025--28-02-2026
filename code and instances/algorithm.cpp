
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "functions.h"

using namespace std;

#define min(X,Y) ((X)<(Y)?(X):(Y))
#define max(X,Y) ((X)>(Y)?(X):(Y))


// Input data:
string data_darp = "datafile_4_request.txt";

// Problem-related parameters:
bool transfers = 1;
bool relative_mrt = 1;
double mrt_factor = 1.50;
double speed_factor_pt = 1.00;
double pt_interval = 10.0;

// Algorithm-related parameters:
int iterations = 25000;
int iter_without_impr = 500;
double max_removal_pct = 0.30;
double max_deterior_pct = 0.04;

// Removal-related parameters:
int clustered_removal = 0;
int idarp_adapted = 1;
int without_operator = 0;

// Insertion-related parameters:
int terminals_per_node = 3;
int priority_non_pt = 0;
int priority_pt = 0;
int priority_long_dist = 0;
int priority_short_dist = 0;

// Variables to measure the efficiency of the operators:
int accept_random_removal = 0, accept_worst_removal = 0, accept_related_removal = 0, accept_route_removal = 0, accept_random_order_insertion = 0, accept_greedy_insertion = 0, accept_two_regret_insertion = 0;
int improve_random_removal = 0, improve_worst_removal = 0, improve_related_removal = 0, improve_route_removal = 0, improve_random_order_insertion = 0, improve_greedy_insertion = 0, improve_two_regret_insertion = 0;


// Main function to run the algorithm:
int main(int argc, char *argv[]) {

	std::cout << "START\n";

	srand(time(NULL));
	chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();

	// Input data via command line:
	if (argc > 1) {
		data_darp =				argv[1];
		transfers =				atoi(argv[2]);
		relative_mrt =			atoi(argv[3]);
		mrt_factor =			atof(argv[4]);
		speed_factor_pt =		atof(argv[5]);
		pt_interval =			atof(argv[6]);
		iterations =			atoi(argv[7]);
		iter_without_impr =		atoi(argv[8]);
		max_removal_pct =		atof(argv[9]);
		max_deterior_pct =		atof(argv[10]);
		clustered_removal =		atoi(argv[11]);
		idarp_adapted =			atoi(argv[12]);
		without_operator =		atoi(argv[13]);
		terminals_per_node =	atoi(argv[14]);
		priority_non_pt =		atoi(argv[15]);
		priority_pt =			atoi(argv[16]);
		priority_long_dist =	atoi(argv[17]);
		priority_short_dist =	atoi(argv[18]);
	}

	// Create and initialize problem:
	struct problem p;
	build_problem(p, data_darp);
	std::cout << "Problem Built\n";
	
	// Create and initialize solutions:
	struct solution s_curr;
	struct solution s_best;
	struct solution s_over;
	initialize_solution(p, s_curr);
	initialize_solution(p, s_best);
	initialize_solution(p, s_over);
	std::cout << "Solution intialised \n";
	std::cout << "Problem with " << p.n_requests << " requests, " << p.n_vehicles << " vehicles, " << p.n_terminals << " terminals \n";
	std::cout << "Nodes: " << p.n_nodes << ", References: " << p.n_references << "\n";

	// Construct an initial solution:
	random_order_insertion(p, s_curr);
	while (check_completeness(p, s_curr) == false) {
		random_removal(p, s_curr, max_removal_pct, 0);
		random_order_insertion(p, s_curr);
		std::cout << s_curr.total_distance << "\n";
	}

	std::cout << "Initial solution is " << s_over.total_distance << "\n"; 

	// Update the best and overall best solution with the initial solution:
	update_solution(p, s_best, s_curr, true);
	update_solution(p, s_over, s_curr, true);
	printf("Best solution %.2f found in initialization \n", s_over.total_distance);
	
	int iter_last_impr = 0;

	// Perform a prefixed number of iterations of the algorithm:
	for (int iteration = 0; iteration < iterations; iteration++) {
		
		// Print intermediate results:
		if (iteration % 5000 == 0) {

			chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();
			chrono::duration<double> difference = end_time - start_time;
			double computation_time = difference.count();

			write_output_file(p, s_over, iteration, computation_time, true);
		}

		// Apply a diversification to s_best if needed:
		if ((iteration - iter_last_impr) == iter_without_impr) {

			do {
				random_removal(p, s_best, max_removal_pct, 0);
				random_order_insertion(p, s_best);
			} while (check_completeness(p, s_best) == false);
			relocate(p, s_best);
			exchange_natural_sequences(p, s_best);

			iter_last_impr = iteration;
			update_solution(p, s_curr, s_best, true);
			printf("Diversification to %.2f in iteration %d \n", s_best.total_distance, iteration + 1);
		}

		// Determine removal percentage and maximum deterioration:
		double removal_pct = 0.01 + ((double)rand() / (RAND_MAX)) * (max_removal_pct - 0.01);
		double max_objective_value = (1 + max_deterior_pct) * s_best.total_distance;
		
		// Perform LNS removal operators on the current solution:
		int d;
		switch (rand() % (4-without_operator)) {
			case 0: random_removal(p, s_curr, removal_pct, clustered_removal); d = 0; break;
			case 1: worst_removal(p, s_curr, removal_pct, clustered_removal, idarp_adapted); d = 1; break;
			case 2: related_removal(p, s_curr, removal_pct, clustered_removal, idarp_adapted); d = 2; break;
			case 3: route_removal(p, s_curr, removal_pct, clustered_removal); d = 3; break;
		}

		// Perform LNS insertion operators on the current solution:
		int r;
		switch (rand() % 3) {
			case 0: random_order_insertion(p, s_curr, max_objective_value); r = 0; break;
			case 1: greedy_insertion(p, s_curr, max_objective_value); r = 1; break;
			case 2: two_regret_insertion(p, s_curr, max_objective_value); r = 2; break;
		}
		
		// If no complete solution has been obtained, return to the previous or the best solution:
		if (check_completeness(p, s_curr) == false) { update_solution(p, s_curr, s_best, true); }

		else {

			if (d == 0) { accept_random_removal += 1; }
			else if (d == 1) { accept_worst_removal += 1; }
			else if (d == 2) { accept_related_removal += 1; }
			else { accept_route_removal += 1; }

			if (r == 0) { accept_random_order_insertion += 1; }
			else if (r == 1) { accept_greedy_insertion += 1; }
			else { accept_two_regret_insertion += 1; }

			// Perform additional local search on complete solutions:
			relocate(p, s_curr);
			exchange_natural_sequences(p, s_curr);
			
			// Double-check the completeness of the solution before updating (not necessary):
			if (check_completeness(p, s_curr) == false ) { update_solution(p, s_curr, s_best, true); printf("ERROR - Incomplete solution corrected\n"); }

			// Update the best solution if needed:
			else if (s_curr.total_distance < s_best.total_distance - 0.001) {

				iter_last_impr = iteration;
				update_solution(p, s_best, s_curr, true);

				// Double-check the feasibility of the solution (not necessary):
				if (check_feasibility(p, s_best) == false) { printf("ERROR - Infeasible best solution accepted\n"); }

				// Update the overall best solution if needed:
				else if (s_curr.total_distance < s_over.total_distance - 0.001) {

					if (d == 0) { improve_random_removal += 1; }
					else if (d == 1) { improve_worst_removal += 1; }
					else if (d == 2) { improve_related_removal += 1; }
					else { improve_route_removal += 1; }

					if (r == 0) { improve_random_order_insertion += 1; }
					else if (r == 1) { improve_greedy_insertion += 1; }
					else { improve_two_regret_insertion += 1; }

					update_solution(p, s_over, s_curr, true);
					printf("Best solution %.2f found after iteration %d \n", s_over.total_distance, iteration + 1);
				}
			}
		}
	}

	// Register total computation time:
	chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();
	chrono::duration<double> difference = end_time - start_time;
	double computation_time = difference.count();

	write_output_file(p, s_over, iterations, computation_time, true);
}


void build_problem(struct problem &p, std::string data_darp) {

	// Read input data:
	read_data_darp(p, data_darp);
	
	// Preprocessing steps:
	compute_matrix_nontightened(p);	
	if (transfers == 1) { compute_matrix_pt(p); }
	compute_relative_mrt(p);
	if (transfers == 1) { compute_nearest_terminals(p); }
	tighten_time_windows(p);
	tighten_matrix(p);
	compute_user_dissimilarities(p);
	compute_node_dissimilarities(p);
}


void read_data_darp(problem &p, std::string data_darp) {

	ifstream infile;
	infile.open(data_darp);
	if (!infile.is_open()) {
		std::cerr << "File not opened\n";
	}


	infile >> p.n_vehicles >> p.n_requests >> p.n_terminals;

	p.n_nodes = 2 * p.n_requests + p.n_vehicles + p.n_terminals;
	p.n_references = 2 * p.n_requests + 2 * p.n_vehicles + 2 * p.n_requests * transfers;

	p.nodes = new node[p.n_nodes];
	p.requests = new request[p.n_requests];
	p.vehicles = new vehicle[p.n_vehicles];

	int node, mobile, wheelchair;
	double x_coord, y_coord, lower_tw, upper_tw, service_dur;

	while (infile >> node >> x_coord >> y_coord >> lower_tw >> upper_tw >> mobile >> wheelchair >> service_dur) {
		
		// Data pickup nodes and requests:
		if (node < p.n_requests) {
			p.nodes[node] = { x_coord, y_coord, lower_tw, upper_tw, service_dur, 0 };

			int request = node;
			p.requests[request] = { mobile, wheelchair, 0.0 };
		}

		// Data delivery nodes:
		else if (node < 2 * p.n_requests) {
			p.nodes[node] = { x_coord, y_coord, lower_tw, upper_tw, service_dur, 1 };
		}

		// Data depot nodes and vehicles:
		else if (node < 2 * p.n_requests + p.n_vehicles) {
			p.nodes[node] = { x_coord, y_coord, lower_tw, upper_tw, service_dur, 2 }; 

			int vehicle = node - 2 * p.n_requests;
			p.vehicles[vehicle] = { mobile, wheelchair, 720.0 };
		}

		// Data transfer nodes:
		else if (transfers == 1) {
			p.nodes[node] = { x_coord, y_coord, lower_tw, upper_tw, service_dur, 3 };
		}
	}

	infile.close();
}


void compute_matrix_nontightened(struct problem &p) {

	p.matrix_nontightened = new double[p.n_nodes * p.n_nodes];

	for (int node1 = 0; node1 < p.n_nodes; node1++) {

		p.matrix_nontightened[node1 * p.n_nodes + node1] = 0.0;

		for (int node2 = 0; node2 < node1; node2++) {
			p.matrix_nontightened[node1 * p.n_nodes + node2] = sqrt(pow(p.nodes[node1].x_coord - p.nodes[node2].x_coord, 2) + pow(p.nodes[node1].y_coord - p.nodes[node2].y_coord, 2));
			p.matrix_nontightened[node2 * p.n_nodes + node1] = p.matrix_nontightened[node1 * p.n_nodes + node2];
		}
	}
}


void compute_matrix_pt(struct problem &p) {

	p.matrix_pt = new double[p.n_terminals * p.n_terminals];

	for (int terminal1 = 0; terminal1 < p.n_terminals; terminal1++) {
		for (int terminal2 = 0; terminal2 <= terminal1; terminal2++) {

			int node_terminal1 = 2 * p.n_requests + p.n_vehicles + terminal1;
			int node_terminal2 = 2 * p.n_requests + p.n_vehicles + terminal2;

			if (terminal1 == terminal2) {
				p.matrix_pt[terminal1 * p.n_terminals + terminal2] = 0;
			}

			// Exclude trips which involve more than one change between PT buses (instances of Molenbruch et al. 2018):
			else if ((p.nodes[node_terminal1].x_coord == -1 * p.nodes[node_terminal2].x_coord && abs(p.nodes[node_terminal1].y_coord) != 10 && abs(p.nodes[node_terminal2].y_coord) != 10) ||
				(p.nodes[node_terminal1].y_coord == -1 * p.nodes[node_terminal2].y_coord && abs(p.nodes[node_terminal1].x_coord) != 10 && abs(p.nodes[node_terminal2].x_coord) != 10)) {
				p.matrix_pt[terminal1 * p.n_terminals + terminal2] = 10000;

			}

			// Compute the travel time between two PT terminals:
			else {
				p.matrix_pt[terminal1 * p.n_terminals + terminal2] = abs(p.nodes[node_terminal1].x_coord - p.nodes[node_terminal2].x_coord) +
					abs(p.nodes[node_terminal1].y_coord - p.nodes[node_terminal2].y_coord);

				// Multiply the travel time by the speed factor of the PT:
				p.matrix_pt[terminal1 * p.n_terminals + terminal2] = p.matrix_pt[terminal1 * p.n_terminals + terminal2] * speed_factor_pt;
			}

			p.matrix_pt[terminal2 * p.n_terminals + terminal1] = p.matrix_pt[terminal1 * p.n_terminals + terminal2];
		}
	}
}


void compute_nearest_terminals(struct problem &p) {

	for (int r = 0; r < p.n_requests; r++) {

		int pickup_node = r;
		int delivery_node = r + p.n_requests;

		double distance_t1[5];
		double distance_t2[5];

		fill(distance_t1, distance_t1 + 5, 10000.0);
		fill(distance_t2, distance_t2 + 5, 10000.0);
		fill(p.requests[r].nearest_t1, p.requests[r].nearest_t1 + 5, -1);
		fill(p.requests[r].nearest_t2, p.requests[r].nearest_t2 + 5, -1);

		// Determine whether the request has a short or long travel distance:
		if (p.matrix_nontightened[pickup_node * p.n_nodes + delivery_node] < 10) {
			p.requests[r].long_dist = false;
		}
		else {
			p.requests[r].long_dist = true;
		}

		// Determine whether the request is eligible for public transport (only ambulant users and sufficient travel distance):
		if (p.requests[r].users_wheelchair > 0) {
			p.requests[r].eligible_pt = false;
		}
		else {
			p.requests[r].eligible_pt = true;

			// Compute the distance between each terminal and the request's pickup/delivery:
			for (int terminal = 0; terminal < p.n_terminals; terminal++) {

				int node_terminal = 2 * p.n_requests + p.n_vehicles + terminal;

				if (p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal] < distance_t1[0]) {
					distance_t1[4] = distance_t1[3]; distance_t1[3] = distance_t1[2]; distance_t1[2] = distance_t1[1]; distance_t1[1] = distance_t1[0];
					p.requests[r].nearest_t1[4] = p.requests[r].nearest_t1[3]; p.requests[r].nearest_t1[3] = p.requests[r].nearest_t1[2]; p.requests[r].nearest_t1[2] = p.requests[r].nearest_t1[1]; p.requests[r].nearest_t1[1] = p.requests[r].nearest_t1[0];
					distance_t1[0] = p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal];
					p.requests[r].nearest_t1[0] = terminal;
				}
				else if (p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal] < distance_t1[1]) {
					distance_t1[4] = distance_t1[3]; distance_t1[3] = distance_t1[2]; distance_t1[2] = distance_t1[1];
					p.requests[r].nearest_t1[4] = p.requests[r].nearest_t1[3]; p.requests[r].nearest_t1[3] = p.requests[r].nearest_t1[2]; p.requests[r].nearest_t1[2] = p.requests[r].nearest_t1[1];
					distance_t1[1] = p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal];
					p.requests[r].nearest_t1[1] = terminal;
				}
				else if (p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal] < distance_t1[2]) {
					distance_t1[4] = distance_t1[3]; distance_t1[3] = distance_t1[2];
					p.requests[r].nearest_t1[4] = p.requests[r].nearest_t1[3]; p.requests[r].nearest_t1[3] = p.requests[r].nearest_t1[2];
					distance_t1[2] = p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal];
					p.requests[r].nearest_t1[2] = terminal;
				}
				else if (p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal] < distance_t1[3]) {
					distance_t1[4] = distance_t1[3];
					p.requests[r].nearest_t1[4] = p.requests[r].nearest_t1[3];
					distance_t1[3] = p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal];
					p.requests[r].nearest_t1[3] = terminal;
				}
				else if (p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal] < distance_t1[4]) {
					distance_t1[4] = p.matrix_nontightened[pickup_node * p.n_nodes + node_terminal];
					p.requests[r].nearest_t1[4] = terminal;
				}

				if (p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node] < distance_t2[0]) {
					distance_t2[4] = distance_t2[3]; distance_t2[3] = distance_t2[2]; distance_t2[2] = distance_t2[1]; distance_t2[1] = distance_t2[0];
					p.requests[r].nearest_t2[4] = p.requests[r].nearest_t2[3]; p.requests[r].nearest_t2[3] = p.requests[r].nearest_t2[2]; p.requests[r].nearest_t2[2] = p.requests[r].nearest_t2[1]; p.requests[r].nearest_t2[1] = p.requests[r].nearest_t2[0];
					distance_t2[0] = p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node];
					p.requests[r].nearest_t2[0] = terminal;
				}
				else if (p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node] < distance_t2[1]) {
					distance_t2[4] = distance_t2[3]; distance_t2[3] = distance_t2[2]; distance_t2[2] = distance_t2[1];
					p.requests[r].nearest_t2[4] = p.requests[r].nearest_t2[3]; p.requests[r].nearest_t2[3] = p.requests[r].nearest_t2[2]; p.requests[r].nearest_t2[2] = p.requests[r].nearest_t2[1];
					distance_t2[1] = p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node];
					p.requests[r].nearest_t2[1] = terminal;
				}
				else if (p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node] < distance_t2[2]) {
					distance_t2[4] = distance_t2[3]; distance_t2[3] = distance_t2[2];
					p.requests[r].nearest_t2[4] = p.requests[r].nearest_t2[3]; p.requests[r].nearest_t2[3] = p.requests[r].nearest_t2[2];
					distance_t2[2] = p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node];
					p.requests[r].nearest_t2[2] = terminal;
				}
				else if (p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node] < distance_t2[3]) {
					distance_t2[4] = distance_t2[3];
					p.requests[r].nearest_t2[4] = p.requests[r].nearest_t2[3];
					distance_t2[3] = p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node];
					p.requests[r].nearest_t2[3] = terminal;
				}
				else if (p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node] < distance_t2[4]) {
					distance_t2[4] = p.matrix_nontightened[node_terminal * p.n_nodes + delivery_node];
					p.requests[r].nearest_t2[4] = terminal;
				}
			}
		}
	}
}


void compute_relative_mrt(struct problem &p) {

	if (relative_mrt == 1) {
		for (int request = 0; request < p.n_requests; request++) {
			p.requests[request].max_ride_time = p.matrix_nontightened[request * p.n_nodes + (request + p.n_requests)] * mrt_factor;
		}
	}
}


void tighten_time_windows(struct problem &p) {

	// Tighten the time windows of the user nodes:
	for (int request = 0; request < p.n_requests; request++) {

		int pickup_node = request;
		int delivery_node = request + p.n_requests;

		// Tighten time windows of non-critical pickup nodes:
		if (p.nodes[pickup_node].lower_tw == 0.0 && p.nodes[pickup_node].upper_tw == 1440.0) {
			p.nodes[pickup_node].lower_tw = max(0.0, p.nodes[delivery_node].lower_tw - p.requests[request].max_ride_time - p.nodes[pickup_node].service_dur);
			p.nodes[pickup_node].upper_tw = min(1440.0, p.nodes[delivery_node].upper_tw - p.matrix_nontightened[pickup_node * p.n_nodes + delivery_node] - p.nodes[pickup_node].service_dur);
		}

		// Tighten time windows of non-critical delivery nodes:
		else if (p.nodes[delivery_node].lower_tw == 0.0 && p.nodes[delivery_node].upper_tw == 1440.0) {
			p.nodes[delivery_node].lower_tw = max(0.0, p.nodes[pickup_node].lower_tw + p.nodes[pickup_node].service_dur + p.matrix_nontightened[pickup_node * p.n_nodes + delivery_node]);
			p.nodes[delivery_node].upper_tw = min(1440.0, p.nodes[pickup_node].upper_tw + p.nodes[pickup_node].service_dur + p.requests[request].max_ride_time);
		}
	}

	// Tighten the time windows of the transfer terminals:
	for (int terminal = 0; terminal < p.n_terminals; terminal++) {

		int transfer_node = 2 * p.n_requests + p.n_vehicles + terminal;

		double lower_tw = 1440.0;
		double upper_tw = 0.0;

		for (int request = 0; request < p.n_requests; request++) {

			int pickup_node = request;
			int delivery_node = request + p.n_requests;

			// Compute the earliest time by which the transfer node can be reached from the pickup:
			if (p.nodes[pickup_node].lower_tw + p.nodes[pickup_node].service_dur + p.matrix_nontightened[pickup_node * p.n_nodes + transfer_node] < lower_tw) {
				lower_tw = p.nodes[pickup_node].lower_tw + p.nodes[pickup_node].service_dur + p.matrix_nontightened[pickup_node * p.n_nodes + transfer_node];
			}

			// Compute the latest time by which the transfer node should be left to reach the delivery:
			if (p.nodes[delivery_node].upper_tw - p.nodes[transfer_node].service_dur - p.matrix_nontightened[transfer_node * p.n_nodes + delivery_node] > upper_tw) {
				upper_tw = p.nodes[delivery_node].upper_tw - p.nodes[transfer_node].service_dur - p.matrix_nontightened[transfer_node * p.n_nodes + delivery_node];
			}
		}

		// Strengthen the time windows of the transfer nodes if possible:
		if (lower_tw > p.nodes[transfer_node].lower_tw) { p.nodes[transfer_node].lower_tw = lower_tw; }
		if (upper_tw < p.nodes[transfer_node].upper_tw) { p.nodes[transfer_node].upper_tw = upper_tw; }
	}
}


void tighten_matrix(struct problem &p) {

	// Use tightening constraints from Cordeau (2006) adapted to the IDARP
	// Note: most tightening constraints do not apply if speed PT > speed DARP

	p.matrix = new double[p.n_nodes * p.n_nodes];

	for (int node1 = 0; node1 < p.n_nodes; node1++) {
		for (int node2 = 0; node2 < p.n_nodes; node2++) {

			if ((node1 == node2 + p.n_requests && p.nodes[node2].type == 0) ||
				(p.nodes[node1].type == 2 && p.nodes[node2].type == 1) ||
				(p.nodes[node2].type == 2 && p.nodes[node1].type == 0) ||
				(p.nodes[node1].lower_tw + p.nodes[node1].service_dur + p.matrix_nontightened[node1 * p.n_nodes + node2] > p.nodes[node2].upper_tw)) {

				p.matrix[node1 * p.n_nodes + node2] = 1000.0;
			}

			else {
				p.matrix[node1 * p.n_nodes + node2] = p.matrix_nontightened[node1 * p.n_nodes + node2];
			}
		}
	}

	for (int request = 0; request < p.n_requests; request++) {

		int pickup = request;
		int delivery = request + p.n_requests;

		for (int node2 = 0; node2 < p.n_nodes; node2++) {

			if (node2 != pickup && node2 != delivery && (p.nodes[node2].type <= 1 || p.nodes[node2].type == 3) &&
				p.matrix_nontightened[pickup * p.n_nodes + node2] + p.nodes[node2].service_dur + p.matrix_nontightened[node2 * p.n_nodes + delivery] > p.requests[request].max_ride_time) {

				p.matrix[pickup * p.n_nodes + node2] = 1000.0;
				p.matrix[node2 * p.n_nodes + delivery] = 1000.0;
			}
		}
	}

	for (int request1 = 0; request1 < p.n_requests; request1++) {
		for (int request2 = 0; request2 < p.n_requests; request2++) {
			if (request1 != request2) {

				int pickup_node1 = request1;
				int pickup_node2 = request2;
				int delivery_node1 = request1 + p.n_requests;
				int delivery_node2 = request2 + p.n_requests;

				// Path (j, i, n+j, n+i):
				double et_2 = max(p.nodes[pickup_node1].lower_tw, p.nodes[pickup_node2].lower_tw + p.nodes[pickup_node2].service_dur + p.matrix[pickup_node2 * p.n_nodes + pickup_node1]);
				double et_3 = max(p.nodes[delivery_node2].lower_tw, et_2 + p.nodes[pickup_node1].service_dur + p.matrix[pickup_node1 * p.n_nodes + delivery_node2]);
				if (et_3 + p.nodes[delivery_node2].service_dur + p.matrix[delivery_node2 * p.n_nodes + delivery_node1] > p.nodes[delivery_node1].upper_tw) {
					p.matrix[pickup_node1 * p.n_nodes + delivery_node2] = 1000.0;
				}

				// Path (i, n+i, j, n+j):
				et_2 = max(p.nodes[delivery_node1].lower_tw, p.nodes[pickup_node1].lower_tw + p.nodes[pickup_node1].service_dur + p.matrix[pickup_node1 * p.n_nodes + delivery_node1]);
				et_3 = max(p.nodes[pickup_node2].lower_tw, et_2 + p.nodes[delivery_node1].service_dur + p.matrix[delivery_node1 * p.n_nodes + pickup_node2]);
				if (et_3 + p.nodes[pickup_node2].service_dur + p.matrix[pickup_node2 * p.n_nodes + delivery_node2] > p.nodes[delivery_node2].upper_tw) {
					p.matrix[delivery_node1 * p.n_nodes + pickup_node2] = 1000.0;
				}

				// Some reductions are only valid if no transfers can be performed (H�ll et al., 2006):
				if (transfers == 0) {

					// Path (i, j, n+i, n+j):
					et_2 = max(p.nodes[pickup_node2].lower_tw, p.nodes[pickup_node1].lower_tw + p.nodes[pickup_node1].service_dur + p.matrix[pickup_node1 * p.n_nodes + pickup_node2]);
					et_3 = max(p.nodes[delivery_node1].lower_tw, et_2 + p.nodes[pickup_node2].service_dur + p.matrix[pickup_node2 * p.n_nodes + delivery_node1]);
					if (et_3 + p.nodes[delivery_node1].service_dur + p.matrix[delivery_node1 * p.n_nodes + delivery_node2] > p.nodes[delivery_node2].upper_tw) {

						// Path (i, j, n+j, n+i):
						et_3 = max(p.nodes[delivery_node2].lower_tw, et_2 + p.nodes[pickup_node2].service_dur + p.matrix[pickup_node2 * p.n_nodes + delivery_node2]);
						if (et_3 + p.nodes[delivery_node2].service_dur + p.matrix[delivery_node2 * p.n_nodes + delivery_node1] > p.nodes[delivery_node1].upper_tw) {
							p.matrix[pickup_node1 * p.n_nodes + pickup_node2] = 1000.0;
						}

						// Path (j, i, n+i, n+j):
						et_2 = max(p.nodes[pickup_node1].lower_tw, p.nodes[pickup_node2].lower_tw + p.nodes[pickup_node2].service_dur + p.matrix[pickup_node2 * p.n_nodes + pickup_node1]);
						et_3 = max(p.nodes[delivery_node1].lower_tw, et_2 + p.nodes[pickup_node1].service_dur + p.matrix[pickup_node1 * p.n_nodes + delivery_node1]);
						if (et_3 + p.nodes[delivery_node1].service_dur + p.matrix[delivery_node1 * p.n_nodes + delivery_node2] > p.nodes[delivery_node2].upper_tw) {
							p.matrix[delivery_node1 * p.n_nodes + delivery_node2] = 1000.0;
						}
					}
				}
			}
		}
	}
}


void compute_user_dissimilarities(struct problem &p) {

	int n_requests = p.n_requests;
	double largest_dist_dissim = 0.0, largest_time_dissim = 0.0;

	double *user_dissim_dist = new double[n_requests * n_requests];
	double *user_dissim_time = new double[n_requests * n_requests];

	for (int request1 = 0; request1 < n_requests; request1++) {

		user_dissim_dist[request1 * n_requests + request1] = 0.0;
		user_dissim_time[request1 * n_requests + request1] = 0.0;

		int pickup1 = request1;
		int delivery1 = pickup1 + n_requests;

		for (int request2 = 0; request2 < request1; request2++) {

			int pickup2 = request2;
			int delivery2 = pickup2 + n_requests;

			double dist_dissim = p.matrix_nontightened[pickup1 * p.n_nodes + pickup2] + p.matrix_nontightened[delivery1 * p.n_nodes + delivery2];

			user_dissim_dist[request1 * n_requests + request2] = dist_dissim;
			user_dissim_dist[request2 * n_requests + request1] = dist_dissim;

			if (dist_dissim > largest_dist_dissim) { largest_dist_dissim = dist_dissim; }

			double time_dissim = abs(p.nodes[pickup1].lower_tw + p.nodes[pickup1].upper_tw - p.nodes[pickup2].lower_tw - p.nodes[pickup2].upper_tw) +
								 abs(p.nodes[delivery1].lower_tw + p.nodes[delivery1].upper_tw - p.nodes[delivery2].lower_tw - p.nodes[delivery2].upper_tw);

			user_dissim_time[request1 * n_requests + request2] = time_dissim;
			user_dissim_time[request2 * n_requests + request1] = time_dissim;

			if (time_dissim > largest_time_dissim) { largest_time_dissim = time_dissim; }
		}
	}

	p.user_dissim = new double[n_requests * n_requests];

	for (int request1 = 0; request1 < n_requests; request1++) {
		for (int request2 = 0; request2 <= request1; request2++) {

			double dissim = user_dissim_dist[request1 * n_requests + request2] / largest_dist_dissim + user_dissim_time[request1 * n_requests + request2] / largest_time_dissim;

			p.user_dissim[request1 * n_requests + request2] = dissim;
			p.user_dissim[request2 * n_requests + request1] = dissim;
		}
	}

	delete[] user_dissim_dist;
	delete[] user_dissim_time;
}


void compute_node_dissimilarities(struct problem &p) {

	int n_requests = p.n_requests;
	int n_user_locations = 2 * p.n_requests;
	double largest_dist_dissim = 0.0, largest_time_dissim = 0.0;

	double *node_dissim_dist = new double[n_user_locations * n_user_locations];
	double *node_dissim_time = new double[n_user_locations * n_user_locations];

	for (int node1 = 0; node1 < n_user_locations; node1++) {

		node_dissim_dist[node1 * n_user_locations + node1] = 0.0;
		node_dissim_time[node1 * n_user_locations + node1] = 0.0;

		for (int node2 = 0; node2 < node1; node2++) {
			
			double dist_dissim = p.matrix_nontightened[node1 * p.n_nodes + node2];

			node_dissim_dist[node1 * n_user_locations + node2] = dist_dissim;
			node_dissim_dist[node2 * n_user_locations + node1] = dist_dissim;

			if (dist_dissim > largest_dist_dissim) { largest_dist_dissim = dist_dissim; }

			double time_dissim = abs(p.nodes[node1].lower_tw + p.nodes[node1].upper_tw - p.nodes[node2].lower_tw - p.nodes[node2].upper_tw);
				
			node_dissim_time[node1 * n_user_locations + node2] = time_dissim;
			node_dissim_time[node2 * n_user_locations + node1] = time_dissim;

			if (time_dissim > largest_time_dissim) { largest_time_dissim = time_dissim; }
		}
	}

	p.node_dissim = new double[n_requests * n_requests];

	for (int request1 = 0; request1 < n_requests; request1++) {
		for (int request2 = 0; request2 <= request1; request2++) {

			int pickup1 = request1, delivery1 = request1 + n_requests;
			int pickup2 = request2, delivery2 = request2 + n_requests;

			double dissim = min(
				min(node_dissim_dist[pickup1 * n_user_locations + pickup2] / largest_dist_dissim + node_dissim_time[pickup1 * n_user_locations + pickup2] / largest_time_dissim,
				    node_dissim_dist[pickup1 * n_user_locations + delivery2] / largest_dist_dissim + node_dissim_time[pickup1 * n_user_locations + delivery2] / largest_time_dissim),
				min(node_dissim_dist[delivery1 * n_user_locations + delivery2] / largest_dist_dissim + node_dissim_time[delivery1 * n_user_locations + delivery2] / largest_time_dissim,
				    node_dissim_dist[delivery1 * n_user_locations + pickup2] / largest_dist_dissim + node_dissim_time[delivery1 * n_user_locations + pickup2] / largest_time_dissim));

			p.node_dissim[request1 * n_requests + request2] = dissim;
			p.node_dissim[request2 * n_requests + request1] = dissim;
		}
	}

	delete[] node_dissim_dist;
	delete[] node_dissim_time;
}


void initialize_solution(struct problem &p, struct solution &s) {
	
	s.total_distance = 0.0;
	s.predecessor = new int[p.n_references];
	s.successor = new int[p.n_references];

	s.map_to_node = new int[p.n_references];
	s.map_to_request = new int[p.n_references];
	s.map_to_vehicle = new int[p.n_references];

	s.load_mobile = new int[p.n_references];
	s.load_wheelchair = new int[p.n_references];

	s.service_start = new double[p.n_references];
	s.earliest_time = new double[p.n_references];
	s.latest_time = new double[p.n_references];

	s.best_insertion_cost = new double[p.n_requests];
	s.second_best_insertion_cost = new double[p.n_requests];
	s.best_vehicle1 = new int[p.n_references];
	s.best_vehicle2 = new int[p.n_references];
	s.second_best_vehicle1 = new int[p.n_references];
	s.second_best_vehicle2 = new int[p.n_references];
	s.best_predecessor = new int[p.n_references];
	s.best_successor = new int[p.n_references];
	s.best_terminal1 = new int[p.n_requests];
	s.best_terminal2 = new int[p.n_requests];

	fill(s.predecessor, s.predecessor + p.n_references, -1);
	fill(s.successor, s.successor + p.n_references, -1);

	fill(s.map_to_node, s.map_to_node + p.n_references, -1);
	fill(s.map_to_request, s.map_to_request + p.n_references, -1);
	fill(s.map_to_vehicle, s.map_to_vehicle + p.n_references, -1);

	fill(s.load_mobile, s.load_mobile + p.n_references, 0);
	fill(s.load_wheelchair, s.load_wheelchair + p.n_references, 0);

	fill(s.service_start, s.service_start + p.n_references, 0.0);
	fill(s.earliest_time, s.earliest_time + p.n_references, 0.0);
	fill(s.latest_time, s.latest_time + p.n_references, 1440.0);

	fill(s.best_insertion_cost, s.best_insertion_cost + p.n_requests, 1000.0);
	fill(s.second_best_insertion_cost, s.second_best_insertion_cost + p.n_requests, 1000.0);
	fill(s.best_vehicle1, s.best_vehicle1 + p.n_requests, -1);
	fill(s.best_vehicle2, s.best_vehicle2 + p.n_requests, -1);
	fill(s.second_best_vehicle1, s.second_best_vehicle1 + p.n_requests, -1);
	fill(s.second_best_vehicle2, s.second_best_vehicle2 + p.n_requests, -1);
	fill(s.best_predecessor, s.best_predecessor + p.n_references, -1);
	fill(s.best_successor, s.best_successor + p.n_references, -1);
	fill(s.best_terminal1, s.best_terminal1 + p.n_requests, -1);
	fill(s.best_terminal2, s.best_terminal2 + p.n_requests, -1);

	// Link the references to their corresponding node, request and vehicle (if not variable):
	for (int reference = 0; reference < p.n_references; reference++) {

		// References to pickup node:
		if (reference < p.n_requests) {
			s.map_to_node[reference] = reference;
			s.map_to_request[reference] = reference;
		}

		// References to delivery node:
		else if (reference < 2 * p.n_requests) {
			s.map_to_node[reference] = reference;
			s.map_to_request[reference] = reference - p.n_requests;
		}

		// References to start depot:
		else if (reference < 2 * p.n_requests + p.n_vehicles) {
			s.map_to_node[reference] = reference;
			s.map_to_vehicle[reference] = reference - 2 * p.n_requests;
		}

		// References to end depot:
		else if (reference < 2 * p.n_requests + 2 * p.n_vehicles) {
			s.map_to_node[reference] = reference - p.n_vehicles;
			s.map_to_vehicle[reference] = reference - 2 * p.n_requests - p.n_vehicles;
		}

		// References to first transfer terminal (split up per request):
		else if (reference < 3 * p.n_requests + 2 * p.n_vehicles) {
			s.map_to_request[reference] = reference - 2 * p.n_requests - 2 * p.n_vehicles;
		}

		// References to second transfer terminal (split up per request):
		else {
			s.map_to_request[reference] = reference - 3 * p.n_requests - 2 * p.n_vehicles;
		}
	}

	// Insert the start and end depots into the routes:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int ref_begin_depot = 2 * p.n_requests + vehicle;
		int ref_end_depot = 2 * p.n_requests + p.n_vehicles + vehicle;

		s.successor[ref_begin_depot] = ref_end_depot;
		s.successor[ref_end_depot] = ref_begin_depot;
		s.predecessor[ref_end_depot] = ref_begin_depot;
		s.predecessor[ref_begin_depot] = ref_end_depot;
	}
}


void update_solution(struct problem &p, struct solution &s1, struct solution &s2, bool create_schedule) {

	s1.total_distance = s2.total_distance;

	for (int reference = 0; reference < p.n_references; reference++) {

		s1.predecessor[reference] = s2.predecessor[reference];
		s1.successor[reference] = s2.successor[reference];

		s1.map_to_node[reference] = s2.map_to_node[reference];
		s1.map_to_request[reference] = s2.map_to_request[reference];
		s1.map_to_vehicle[reference] = s2.map_to_vehicle[reference];

		s1.load_mobile[reference] = s2.load_mobile[reference];
		s1.load_wheelchair[reference] = s2.load_wheelchair[reference];

		s1.service_start[reference] = s2.service_start[reference];
		s1.earliest_time[reference] = s2.earliest_time[reference];
		s1.latest_time[reference] = s2.latest_time[reference];
	}

	// Update time schedule:
	if (create_schedule == true) {
		if (transfers == 1) { check_service_start_BellmanFord(p, s1, true); }
		if (transfers == 0) { for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) { check_service_start_Tang(p, s1, vehicle, true); } }
	}
}


void random_removal(struct problem &p, struct solution &s, double removal_percentage, bool clustered_removal) {

	int number_removed = 0;

	while (number_removed < p.n_requests * removal_percentage) {

		// Remove random request:
		int request_id = rand() % p.n_requests;
		
		if (s.map_to_vehicle[request_id] != -1) {

			// Remove requests using the same terminals around the same time (function returns increase of number removed);
			if (transfers == 1 && clustered_removal == 1) { number_removed += clear_terminal(p, s, request_id); }

			remove_request(p, s, request_id);
			number_removed += 1;
		}
	}

	// Update the load, earliest time and latest time array for each route:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int begin_depot = 2 * p.n_requests + vehicle;
		int end_depot = begin_depot + p.n_vehicles;

		update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
		update_earliest_time(p, s, begin_depot);
		update_latest_time(p, s, end_depot);
	}
}


void worst_removal(struct problem &p, struct solution &s, double removal_percentage, bool clustered_removal, bool idarp_adapted) {

	double *removal_benefit = new double[p.n_requests];

	int number_removed = 0;

	// First compute the cost benefit of removing a request from the solution (large value = large benefit):
	for (int request = 0; request < p.n_requests; request++) {

		int ref_pickup = request; int node_pickup = s.map_to_node[ref_pickup];
		int ref_delivery = request + p.n_requests; int node_delivery = s.map_to_node[ref_delivery];

		double removal_benefit_pickup = 0;
		double removal_benefit_delivery = 0;

		// Situation 1: the request is partly served by public transport:
		if (s.map_to_vehicle[ref_pickup] != s.map_to_vehicle[ref_delivery]) {

			int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request; int node_terminal1 = s.map_to_node[ref_terminal1];
			int ref_terminal2 = 3 * p.n_requests + 2 * p.n_vehicles + request; int node_terminal2 = s.map_to_node[ref_terminal2];

			if (ref_terminal1 == s.successor[ref_pickup]) {

				int ref_predecessor_pickup = s.predecessor[ref_pickup];
				int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
				int ref_successor_terminal1 = s.successor[ref_terminal1];
				int node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

				removal_benefit_pickup = p.matrix[node_predecessor_pickup * p.n_nodes + node_pickup] +
										 p.matrix[node_pickup * p.n_nodes + node_terminal1] +
										 p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1] -
										 p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_terminal1];
			}

			else {

				int ref_predecessor_pickup = s.predecessor[ref_pickup];
				int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
				int ref_successor_pickup = s.successor[ref_pickup];
				int node_successor_pickup = s.map_to_node[ref_successor_pickup];
				int ref_predecessor_terminal1 = s.predecessor[ref_terminal1];
				int node_predecessor_terminal1 = s.map_to_node[ref_predecessor_terminal1];
				int ref_successor_terminal1 = s.successor[ref_terminal1];
				int node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

				removal_benefit_pickup = p.matrix[node_predecessor_pickup * p.n_nodes + node_pickup] +
										 p.matrix[node_pickup * p.n_nodes + node_successor_pickup] +
										 p.matrix[node_predecessor_terminal1 * p.n_nodes + node_terminal1] +
										 p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1] -
										 p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup] -
										 p.matrix[node_predecessor_terminal1 * p.n_nodes + node_successor_terminal1];
			}

			if (ref_delivery == s.successor[ref_terminal2]) {

				int ref_predecessor_terminal2 = s.predecessor[ref_terminal2];
				int node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];
				int ref_successor_delivery = s.successor[ref_delivery];
				int node_successor_delivery = s.map_to_node[ref_successor_delivery];

				removal_benefit_delivery = p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] +
										   p.matrix[node_terminal2 * p.n_nodes + node_delivery] +
										   p.matrix[node_delivery * p.n_nodes + node_successor_delivery] -
										   p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_delivery];
			}

			else {

				int ref_predecessor_terminal2 = s.predecessor[ref_terminal2];
				int node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];
				int ref_successor_terminal2 = s.successor[ref_terminal2];
				int node_successor_terminal2 = s.map_to_node[ref_successor_terminal2];
				int ref_predecessor_delivery = s.predecessor[ref_delivery];
				int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
				int ref_successor_delivery = s.successor[ref_delivery];
				int node_successor_delivery = s.map_to_node[ref_successor_delivery];

				removal_benefit_delivery = p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] +
										   p.matrix[node_terminal2 * p.n_nodes + node_successor_terminal2] +
										   p.matrix[node_predecessor_delivery * p.n_nodes + node_delivery] +
										   p.matrix[node_delivery * p.n_nodes + node_successor_delivery] -
										   p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_terminal2] -
										   p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];
			}
		}

		// Situation 2: the request is not served by public transport:
		else {

			if (ref_delivery == s.successor[ref_pickup]) {

				int ref_predecessor_pickup = s.predecessor[ref_pickup];
				int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
				int ref_successor_delivery = s.successor[ref_delivery];
				int node_successor_delivery = s.map_to_node[ref_successor_delivery];

				double removal_benefit = p.matrix[node_predecessor_pickup * p.n_nodes + node_pickup] +
										 p.matrix[node_pickup * p.n_nodes + node_delivery] +
										 p.matrix[node_delivery * p.n_nodes + node_successor_delivery] -
										 p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_delivery];

				removal_benefit_pickup = removal_benefit / 2;
				removal_benefit_delivery = removal_benefit / 2;
			}

			else {

				int ref_predecessor_pickup = s.predecessor[ref_pickup];
				int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
				int ref_successor_pickup = s.successor[ref_pickup];
				int node_successor_pickup = s.map_to_node[ref_successor_pickup];
				int ref_predecessor_delivery = s.predecessor[ref_delivery];
				int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
				int ref_successor_delivery = s.successor[ref_delivery];
				int node_successor_delivery = s.map_to_node[ref_successor_delivery];

				removal_benefit_pickup = p.matrix[node_predecessor_pickup * p.n_nodes + node_pickup] +
										 p.matrix[node_pickup * p.n_nodes + node_successor_pickup] -
										 p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup];

				removal_benefit_delivery = p.matrix[node_predecessor_delivery * p.n_nodes + node_delivery] +
										   p.matrix[node_delivery * p.n_nodes + node_successor_delivery] -
										   p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];
			}
		}

		if (idarp_adapted == 1) { removal_benefit[request] = max(removal_benefit_pickup, removal_benefit_delivery); }
		else { removal_benefit[request] = removal_benefit_pickup + removal_benefit_delivery; }

		// Apply random noise:
		double noise = 0.7 + (double)rand() / RAND_MAX * (1.3 - 0.7);
		removal_benefit[request] = noise * removal_benefit[request];
	}

	// Continue removing requests until the required number is reached:
	while (number_removed < p.n_requests * removal_percentage) {

		double largest_benefit = 0.001;
		int best_request = -1;

		// Select the request for which the removal benefit is largest:
		for (int request = 0; request < p.n_requests; request++) {
			if (removal_benefit[request] > largest_benefit) {
				largest_benefit = removal_benefit[request];
				best_request = request;
			}
		}

		// Remove requests using the same terminals around the same time (adapted version of "clear terminal" function):
		if (transfers == 1 && clustered_removal == 1) { 
			
			int ref_request1_terminal1 = best_request + 2 * p.n_requests + 2 * p.n_vehicles;
			int ref_request1_terminal2 = ref_request1_terminal1 + p.n_requests;

			if (s.map_to_vehicle[ref_request1_terminal1] != -1) {

				for (int request2_id = 0; request2_id < p.n_requests; request2_id++) {

					int ref_request2_terminal1 = request2_id + 2 * p.n_requests + 2 * p.n_vehicles;
					int ref_request2_terminal2 = ref_request2_terminal1 + p.n_requests;

					if (request2_id != best_request && s.map_to_vehicle[ref_request2_terminal1] != -1 &&
						((s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal1] + 30) ||
						(s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal1] + 30) ||
						(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal2] + 30) ||
						(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal2] + 30))) {
						
						remove_request(p, s, request2_id);
						number_removed += 1;
						removal_benefit[request2_id] = 0.0;
					}
				}
			}
		}

		remove_request(p, s, best_request);
		number_removed += 1;
		removal_benefit[best_request] = 0.0;
	}

	// Update the load, earliest time and latest time array for each route:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int begin_depot = 2 * p.n_requests + vehicle;
		int end_depot = begin_depot + p.n_vehicles;

		update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
		update_earliest_time(p, s, begin_depot);
		update_latest_time(p, s, end_depot);
	}

	delete[] removal_benefit;
}


void related_removal(struct problem &p, struct solution &s, double removal_percentage, bool clustered_removal, bool idarp_adapted) {

	int number_removed = 0;
	double *dissimilarity = new double[p.n_requests];

	// Randomly select a first candidate request to be removed:
	int first_request = rand() % (p.n_requests);

	remove_request(p, s, first_request);
	number_removed += 1;
	dissimilarity[first_request] = 1000.0;

	// Compute the dissimilarity of each other request to this first request (small value = highly similar):
	for (int request = 0; request < p.n_requests; request++) {
		if (request != first_request) {

			if (idarp_adapted == 1) { dissimilarity[request] = p.node_dissim[first_request * p.n_requests + request]; }
			else { dissimilarity[request] = p.user_dissim[first_request * p.n_requests + request]; }

			// Apply random noise to avoid that the same requests are often removed:
			double noise = 0.7 + (double)rand() / RAND_MAX * (1.3 - 0.7);
			dissimilarity[request] = noise * dissimilarity[request];
		}
	}

	// Continue removing requests until the required number is reached:
	while (number_removed < p.n_requests * removal_percentage) {

		double best_dissim = 999.9;
		int best_request = -1;

		// Remove the request that is most related to the first request:
		for (int request = 0; request < p.n_requests; request++) {
			if (dissimilarity[request] < best_dissim) {
				best_dissim = dissimilarity[request];
				best_request = request;
			}
		}

		// Remove requests using the same terminals around the same time (adapted version of "clear terminal" function):
		if (transfers == 1 && clustered_removal == 1) {

			int ref_request1_terminal1 = best_request + 2 * p.n_requests + 2 * p.n_vehicles;
			int ref_request1_terminal2 = ref_request1_terminal1 + p.n_requests;

			if (s.map_to_vehicle[ref_request1_terminal1] != -1) {

				for (int request2_id = 0; request2_id < p.n_requests; request2_id++) {

					int ref_request2_terminal1 = request2_id + 2 * p.n_requests + 2 * p.n_vehicles;
					int ref_request2_terminal2 = ref_request2_terminal1 + p.n_requests;

					if (request2_id != best_request && s.map_to_vehicle[ref_request2_terminal1] != -1 &&
						((s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal1] + 30) ||
						(s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal1] + 30) ||
						(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal2] + 30) ||
						(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal2] + 30))) {

						remove_request(p, s, request2_id);
						number_removed += 1;
						dissimilarity[request2_id] = 1000.0;
					}
				}
			}
		}

		remove_request(p, s, best_request);
		number_removed += 1;
		dissimilarity[best_request] = 1000.0;
	}

	// Update the load, earliest time and latest time array for each route:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int begin_depot = 2 * p.n_requests + vehicle;
		int end_depot = begin_depot + p.n_vehicles;

		update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
		update_earliest_time(p, s, begin_depot);
		update_latest_time(p, s, end_depot);
	}

	delete[] dissimilarity;
}


void route_removal(struct problem &p, struct solution &s, double removal_percentage, bool clustered_removal) {

	int number_removed = 0;

	while (number_removed < p.n_requests * removal_percentage) {

		// Select a random vehicle:
		int vehicle = rand() % (p.n_vehicles);
		int ref_begin_depot = 2 * p.n_requests + vehicle;
		int ref_end_depot = ref_begin_depot + p.n_vehicles;

		// Remove requests that are served by this vehicle:
		if (ref_begin_depot != s.predecessor[ref_end_depot]) {

			for (int request_id = 0; request_id < p.n_requests; request_id++) {
				if (s.map_to_vehicle[request_id] == vehicle) {

					// Remove requests using the same terminals around the same time (function returns increase of number removed);
					if (transfers == 1 && clustered_removal == 1) { number_removed += clear_terminal(p, s, request_id); }

					remove_request(p, s, request_id);
					number_removed += 1;
				}
			}
		}
	}

	// Update the load, earliest time and latest time array for each route:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int begin_depot = 2 * p.n_requests + vehicle;
		int end_depot = begin_depot + p.n_vehicles;

		update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
		update_earliest_time(p, s, begin_depot);
		update_latest_time(p, s, end_depot);
	}
}


void clear_insertion_data(struct problem &p, struct solution &s) {

	fill(s.best_insertion_cost, s.best_insertion_cost + p.n_requests, 1000.0);
	fill(s.second_best_insertion_cost, s.second_best_insertion_cost + p.n_requests, 1000.0);
	fill(s.best_vehicle1, s.best_vehicle1 + p.n_requests, -1);
	fill(s.best_vehicle2, s.best_vehicle2 + p.n_requests, -1);
	fill(s.second_best_vehicle1, s.second_best_vehicle1 + p.n_requests, -1);
	fill(s.second_best_vehicle2, s.second_best_vehicle2 + p.n_requests, -1);
	fill(s.best_predecessor, s.best_predecessor + p.n_references, -1);
	fill(s.best_successor, s.best_successor + p.n_references, -1);
	fill(s.best_terminal1, s.best_terminal1 + p.n_requests, -1);
	fill(s.best_terminal2, s.best_terminal2 + p.n_requests, -1);

}


void construct_request_banks(struct problem &p, struct solution &s) {

	s.request_bank_prior = {};
	s.request_bank_non_prior = {};

	bool priorities = 0;
	if (rand() / double(RAND_MAX) < 0.50 && accept_random_order_insertion != 0) { priorities = 1; }

	// Assign uninserted requests to priority or non-priority request bank:
	for (int request = 0; request < p.n_requests; request++) {
		if (s.map_to_vehicle[request] == -1) {

			// Priority requests:
			if (priorities == 1 &&
				((priority_pt == 1 && p.requests[request].eligible_pt == 1) ||
				(priority_non_pt == 1 && p.requests[request].eligible_pt == 0) ||
				(priority_long_dist == 1 && p.requests[request].long_dist == 1) ||
				(priority_short_dist == 1 && p.requests[request].long_dist == 0))) {
				s.request_bank_prior.push_back(request);
			}

			// Non-priority requests:
			else {
				s.request_bank_non_prior.push_back(request);
			}
		}
	}
}


void compute_best_insertion_request(struct problem &p, struct solution &s, int request, double max_objective_value, bool track_second_best) {

	s.best_insertion_cost[request] = min(999.9, max_objective_value - s.total_distance);
	s.second_best_insertion_cost[request] = min(999.9, max_objective_value - s.total_distance);

	// Scenario 1 - Without use of public transport:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {
		compute_best_insertion_without_pt(p, s, request, vehicle, track_second_best);
	}

	// Scenario 2 - With use of public transport:
	if (transfers == 1 && p.requests[request].eligible_pt == true && p.requests[request].long_dist == true) {
		for (int vehicle1 = 0; vehicle1 < p.n_vehicles; vehicle1++) {
			for (int vehicle2 = 0; vehicle2 < p.n_vehicles; vehicle2++) {
				if (vehicle1 != vehicle2) {
					for (int index_t1 = 0; index_t1 < terminals_per_node; index_t1++) {
						int terminal1 = p.requests[request].nearest_t1[index_t1];
						for (int index_t2 = 0; index_t2 < terminals_per_node; index_t2++) {
							int terminal2 = p.requests[request].nearest_t2[index_t2];
							if (terminal1 != terminal2) {
								compute_best_insertion_pt(p, s, request, vehicle1, vehicle2, terminal1, terminal2, track_second_best);
							}
						}
					}
				}
			}
		}
	}
}


void update_best_insertion_request(struct problem &p, struct solution &s, int request, int update_vehicle1, int update_vehicle2, double max_objective_value, bool track_second_best) {

	//////////////////////////////////////////////////////////////////////////////////////////
	// CODE 1: the entire computation of insertion costs is repeated after every insertion: //
	//////////////////////////////////////////////////////////////////////////////////////////

	//s.best_insertion_cost[request] = 1000.0;
	//s.second_best_insertion_cost[request] = 1000.0;
	//s.best_vehicle1[request] = -1;
	//s.best_vehicle2[request] = -1;
	//s.second_best_vehicle1[request] = -1;
	//s.second_best_vehicle2[request] = -1;

	//compute_best_insertion_request(p, s, request, max_objective_value, track_second_best);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// CODE 2: the computation of insertion costs is only repeated for the routes that have changed //
	// (this ignores the fact that an insertion in a route may become infeasible due to a change in //
	//  another route, which will be detected by the feasibility check after the iteration anyway): //
	//////////////////////////////////////////////////////////////////////////////////////////////////

	// Recompute the insertion possibilities if a requests' preferred route has changed:
	if (s.best_vehicle1[request] != -1 &&
		(s.best_vehicle1[request] == update_vehicle1 || s.best_vehicle2[request] == update_vehicle1 ||
		 s.best_vehicle1[request] == update_vehicle2 || (update_vehicle2 != -1 && s.best_vehicle2[request] == update_vehicle2))) {

		// Reset the insertion matrices only for the request under consideration:
		s.best_insertion_cost[request] = 1000.0;
		s.second_best_insertion_cost[request] = 1000.0;
		s.best_vehicle1[request] = -1;
		s.best_vehicle2[request] = -1;
		s.second_best_vehicle1[request] = -1;
		s.second_best_vehicle2[request] = -1;

		compute_best_insertion_request(p, s, request, max_objective_value, track_second_best);
	}

	// If other routes have changed, check for better insertion positions than before in these routes only:
	else if (s.best_vehicle1[request] != -1) {

		compute_best_insertion_without_pt(p, s, request, update_vehicle1, track_second_best);
		if (update_vehicle2 != -1) { compute_best_insertion_without_pt(p, s, request, update_vehicle2, track_second_best); }

		if (transfers == 1 && p.requests[request].eligible_pt == true && p.requests[request].long_dist == true) {
			for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {
				if (vehicle != update_vehicle1) {
					for (int index_t1 = 0; index_t1 < terminals_per_node; index_t1++) {
						int terminal1 = p.requests[request].nearest_t1[index_t1];
						for (int index_t2 = 0; index_t2 < terminals_per_node; index_t2++) {
							int terminal2 = p.requests[request].nearest_t2[index_t2];
							if (terminal1 != terminal2) {
								compute_best_insertion_pt(p, s, request, vehicle, update_vehicle1, terminal1, terminal2, track_second_best);
								compute_best_insertion_pt(p, s, request, update_vehicle1, vehicle, terminal1, terminal2, track_second_best);
							}
						}
					}
				}
				if (vehicle != update_vehicle2 && update_vehicle2 != -1) {
					for (int index_t1 = 0; index_t1 < terminals_per_node; index_t1++) {
						int terminal1 = p.requests[request].nearest_t1[index_t1];
						for (int index_t2 = 0; index_t2 < terminals_per_node; index_t2++) {
							int terminal2 = p.requests[request].nearest_t2[index_t2];
							if (terminal1 != terminal2) {
								compute_best_insertion_pt(p, s, request, vehicle, update_vehicle2, terminal1, terminal2, track_second_best);
								compute_best_insertion_pt(p, s, request, update_vehicle2, vehicle, terminal1, terminal2, track_second_best);
							}
						}
					}
				}
			}
		}
	}
}


void random_order_insertion(struct problem &p, struct solution &s, double max_objective_value) {
	
	clear_insertion_data(p, s);
	construct_request_banks(p, s);

	vector<int> request_bank = {};

	// Randomly order the requests within both priority groups:
	random_shuffle(s.request_bank_prior.begin(), s.request_bank_prior.end());
	random_shuffle(s.request_bank_non_prior.begin(), s.request_bank_non_prior.end());

	// Concatenate both vectors in the right order:
	request_bank.insert(request_bank.end(), s.request_bank_prior.begin(), s.request_bank_prior.end());
	request_bank.insert(request_bank.end(), s.request_bank_non_prior.begin(), s.request_bank_non_prior.end());

	// Select the next request to be inserted:
	for (int index = 0; index < request_bank.size(); index++) {

		int request = request_bank[index];
		compute_best_insertion_request(p, s, request, max_objective_value, false);
		
		// Perform the best feasible insertion (if any):
		if (s.best_vehicle2[request] != -1) { insert_request_with_transfer(p, s, request, s.best_vehicle1[request], s.best_vehicle2[request]); }
		else if (s.best_vehicle1[request] != -1) { insert_request(p, s, request, s.best_vehicle1[request]); }
		else { break; }
	}
}


void greedy_insertion(struct problem &p, struct solution &s, double max_objective_value) {

	clear_insertion_data(p, s);
	construct_request_banks(p, s);

	for (int prior = 1; prior >= 0; prior--) {

		vector<int> request_bank = {};

		// Consider the right request bank based on priorities:
		if (prior == 1) {
			for (int index = 0; index < s.request_bank_prior.size(); index++) {
				request_bank.push_back(s.request_bank_prior[index]);
			}
		}
		else {
			for (int index = 0; index < s.request_bank_non_prior.size(); index++) {
				request_bank.push_back(s.request_bank_non_prior[index]);
			}
		}
		
		// Compute the insertion costs for the requests in the request bank:
		for (int index = 0; index < request_bank.size(); index++) {

			int request = request_bank[index];
			compute_best_insertion_request(p, s, request, max_objective_value, false);
		}

		// Repeat the insertion procedure until no more insertions are possible:
		while (true) {

			int best_request = -1;
			double best_cost = min(999.9, max_objective_value - s.total_distance);

			// Determine the cheapest insertion:
			for (int index = 0; index < request_bank.size(); index++) {

				int request = request_bank[index];
				if (s.map_to_vehicle[request] == -1) {
					if (s.best_insertion_cost[request] < best_cost - 0.001) {

						best_request = request;
						best_cost = s.best_insertion_cost[request];
					}
				}
			}
			
			// Perform the best feasible insertion (if any):
			if (best_request != -1) {
				if (s.best_vehicle2[best_request] != -1) { insert_request_with_transfer(p, s, best_request, s.best_vehicle1[best_request], s.best_vehicle2[best_request]); }
				else { insert_request(p, s, best_request, s.best_vehicle1[best_request]); }
			}
			else { break; }

			// Update the insertion matrix after each insertion:
			for (int index = 0; index < request_bank.size(); index++) {

				int request = request_bank[index];
				if (s.map_to_vehicle[request] == -1) {

					update_best_insertion_request(p, s, request, s.best_vehicle1[best_request], s.best_vehicle2[best_request], max_objective_value, false);

					if (s.best_vehicle1[request] == -1) { return; }
				}
			}
		}
	}
}


void two_regret_insertion(struct problem &p, struct solution &s, double max_objective_value) {

	clear_insertion_data(p, s);
	construct_request_banks(p, s);

	for (int prior = 1; prior >= 0; prior--) {

		vector<int> request_bank = {};

		// Consider the right request bank based on priorities:
		if (prior == 1) {
			for (int index = 0; index < s.request_bank_prior.size(); index++) {
				request_bank.push_back(s.request_bank_prior[index]);
			}
		}
		else {
			for (int index = 0; index < s.request_bank_non_prior.size(); index++) {
				request_bank.push_back(s.request_bank_non_prior[index]);
			}
		}

		// Compute the insertion costs for the requests in the request bank:
		for (int index = 0; index < request_bank.size(); index++) {

			int request = request_bank[index];
			compute_best_insertion_request(p, s, request, max_objective_value, true);
		}

		// Repeat the insertion procedure until no more insertions are possible:
		while (true) {

			int best_request = -1;
			double largest_regret = -0.001;

			// Determine the largest regret over all requests:
			for (int index = 0; index < request_bank.size(); index++) {

				int request = request_bank[index];
				if (s.map_to_vehicle[request] == -1 && s.best_vehicle1[request] != -1) {

					// Update the regret value of each request:
					double regret_value = s.second_best_insertion_cost[request] - s.best_insertion_cost[request];
					if (regret_value > largest_regret) {

						largest_regret = regret_value;
						best_request = request;
					}
				}
			}
			
			// Perform the best feasible insertion (if any):
			if (best_request != -1) {
				if (s.best_vehicle2[best_request] != -1) { insert_request_with_transfer(p, s, best_request, s.best_vehicle1[best_request], s.best_vehicle2[best_request]); }
				else { insert_request(p, s, best_request, s.best_vehicle1[best_request]); }
			}
			else { break; }

			// Update the insertion matrix after each insertion:
			for (int index = 0; index < request_bank.size(); index++) {

				int request = request_bank[index];
				if (s.map_to_vehicle[request] == -1) {

					update_best_insertion_request(p, s, request, s.best_vehicle1[best_request], s.best_vehicle2[best_request], max_objective_value, true);

					if (s.best_vehicle1[request] == -1) { return; }
				}
			}
		}
	}
}


void relocate(struct problem &p, struct solution &s) {

	clear_insertion_data(p, s);

	vector<int> all_requests;

	// Construct an array of all requests:
	for (int request = 0; request < p.n_requests; request++) {
		all_requests.push_back(request);
	}

	// Randomly order the requests:
	random_shuffle(all_requests.begin(), all_requests.end());

	// Select the next request to be removed and reinserted:
	for (int index = 0; index < p.n_requests; index++) {

		int request = all_requests[index];
		int vehicle1 = s.map_to_vehicle[request];
		int vehicle2 = s.map_to_vehicle[request+p.n_requests];

		double max_objective_value = s.total_distance + 0.002;

		remove_request(p, s, request);

		// Update the load, earliest time and latest time array for each route involved:
		int begin_depot = 2 * p.n_requests + vehicle1;
		int end_depot = begin_depot + p.n_vehicles;
		update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
		update_earliest_time(p, s, begin_depot);
		update_latest_time(p, s, end_depot);

		if (vehicle2 != vehicle1) {
			int begin_depot = 2 * p.n_requests + vehicle2;
			int end_depot = begin_depot + p.n_vehicles;
			update_load(p, s, s.successor[begin_depot], s.predecessor[end_depot]);
			update_earliest_time(p, s, begin_depot);
			update_latest_time(p, s, end_depot);
		}
		
		// Compute the best insertion position for this request:
		compute_best_insertion_request(p, s, request, max_objective_value, false);

		// Perform the best feasible insertion (if any):
		if (s.best_vehicle2[request] != -1) { insert_request_with_transfer(p, s, request, s.best_vehicle1[request], s.best_vehicle2[request]); }
		else if (s.best_vehicle1[request] != -1) { insert_request(p, s, request, s.best_vehicle1[request]); }
		else { break; }
	}
}


void exchange_natural_sequences(struct problem &p, struct solution &s) {

	while (true) {

		double best_exchange_cost = 0.0;

		int best_vehicle_1 = -1, best_vehicle_2 = -1;
		int best_predecessor_cut_1a = -1, best_predecessor_cut_2a = -1;
		int best_predecessor_cut_1b = -1, best_predecessor_cut_2b = -1;
		int best_successor_cut_1a = -1, best_successor_cut_2a = -1;
		int best_successor_cut_1b = -1, best_successor_cut_2b = -1;

		// Loop 1 - Select first natural sequence:
		for (int vehicle_1 = 1; vehicle_1 < p.n_vehicles; vehicle_1++) {

			int ref_begin_depot_1 = 2 * p.n_requests + vehicle_1;
			int ref_end_depot_1 = ref_begin_depot_1 + p.n_vehicles;

			int ref_predecessor_cut_1a = ref_begin_depot_1;
			int node_predecessor_cut_1a = s.map_to_node[ref_predecessor_cut_1a];
			int ref_successor_cut_1a = s.successor[ref_predecessor_cut_1a];
			int node_successor_cut_1a = s.map_to_node[ref_successor_cut_1a];

			if (ref_end_depot_1 != s.successor[ref_begin_depot_1]) {

				while (true) {

					if (s.load_mobile[ref_predecessor_cut_1a] == 0 && s.load_wheelchair[ref_predecessor_cut_1a] == 0) {

						int ref_predecessor_cut_1b = ref_successor_cut_1a;
						int node_predecessor_cut_1b = s.map_to_node[ref_predecessor_cut_1b];
						int ref_successor_cut_1b = s.successor[ref_predecessor_cut_1b];
						int node_successor_cut_1b = s.map_to_node[ref_successor_cut_1b];

						while (true) {
							if (s.load_mobile[ref_predecessor_cut_1b] == 0 && s.load_wheelchair[ref_predecessor_cut_1b] == 0) {

								// Loop 2 - Select second natural sequence:
								for (int vehicle_2 = 0; vehicle_2 < vehicle_1; vehicle_2++) {

									int ref_begin_depot_2 = 2 * p.n_requests + vehicle_2;
									int ref_end_depot_2 = ref_begin_depot_2 + p.n_vehicles;

									int ref_predecessor_cut_2a = ref_begin_depot_2;
									int node_predecessor_cut_2a = s.map_to_node[ref_predecessor_cut_2a];
									int ref_successor_cut_2a = s.successor[ref_predecessor_cut_2a];
									int node_successor_cut_2a = s.map_to_node[ref_successor_cut_2a];

									if (ref_end_depot_2 != s.successor[ref_begin_depot_2]) {

										while (true) {
											if (s.load_mobile[ref_predecessor_cut_2a] == 0 && s.load_wheelchair[ref_predecessor_cut_2a] == 0) {

												// Stop moving predecessor_cut_2a as soon as successor_cut_1 cannot be reached in time:
												if (s.earliest_time[ref_predecessor_cut_2a] + p.nodes[node_predecessor_cut_2a].service_dur + p.matrix_nontightened[node_predecessor_cut_2a * p.n_nodes + node_successor_cut_1a] > p.nodes[node_successor_cut_1a].upper_tw) { break; }

												// Check whether the first node of both sequences can be reached in time:
												if (s.earliest_time[ref_predecessor_cut_1a] + p.nodes[node_predecessor_cut_1a].service_dur + p.matrix[node_predecessor_cut_1a * p.n_nodes + node_successor_cut_2a] <= p.nodes[node_successor_cut_2a].upper_tw
													&& s.earliest_time[ref_predecessor_cut_2a] + p.nodes[node_predecessor_cut_2a].service_dur + p.matrix[node_predecessor_cut_2a * p.n_nodes + node_successor_cut_1a] <= p.nodes[node_successor_cut_1a].upper_tw) {

													int ref_predecessor_cut_2b = ref_successor_cut_2a;
													int node_predecessor_cut_2b = s.map_to_node[ref_predecessor_cut_2b];
													int ref_successor_cut_2b = s.successor[ref_predecessor_cut_2b];
													int node_successor_cut_2b = s.map_to_node[ref_successor_cut_2b];

													while (true) {
														if (s.load_mobile[ref_predecessor_cut_2b] == 0 && s.load_wheelchair[ref_predecessor_cut_2b] == 0) {

															// Check whether the last node of both sequences can be left in time:
															if (s.latest_time[ref_successor_cut_1b] - p.matrix[node_predecessor_cut_2b * p.n_nodes + node_successor_cut_1b] - p.nodes[node_predecessor_cut_2b].service_dur >= p.nodes[node_predecessor_cut_2b].lower_tw
																&& s.latest_time[ref_successor_cut_2b] - p.matrix[node_predecessor_cut_1b * p.n_nodes + node_successor_cut_2b] - p.nodes[node_predecessor_cut_1b].service_dur >= p.nodes[node_predecessor_cut_1b].lower_tw) {

																// Check whether an exchange would reduce the cost:
																double exchange_cost = -p.matrix[node_predecessor_cut_1a * p.n_nodes + node_successor_cut_1a] - p.matrix[node_predecessor_cut_1b * p.n_nodes + node_successor_cut_1b]
																	- p.matrix[node_predecessor_cut_2a * p.n_nodes + node_successor_cut_2a] - p.matrix[node_predecessor_cut_2b * p.n_nodes + node_successor_cut_2b]
																	+ p.matrix[node_predecessor_cut_1a * p.n_nodes + node_successor_cut_2a] + p.matrix[node_predecessor_cut_1b * p.n_nodes + node_successor_cut_2b]
																	+ p.matrix[node_predecessor_cut_2a * p.n_nodes + node_successor_cut_1a] + p.matrix[node_predecessor_cut_2b * p.n_nodes + node_successor_cut_1b];

																if (exchange_cost < best_exchange_cost - 0.001) {

																	int ref = ref_successor_cut_1a;
																	while (true) {
																		s.map_to_vehicle[ref] = vehicle_2;
																		ref = s.successor[ref];
																		if (ref == ref_successor_cut_1b) { break; }
																	}

																	ref = ref_successor_cut_2a;
																	while (true) {
																		s.map_to_vehicle[ref] = vehicle_1;
																		ref = s.successor[ref];
																		if (ref == ref_successor_cut_2b) { break; }
																	}

																	// Temporaly perform the exchange:
																	s.predecessor[ref_successor_cut_1a] = ref_predecessor_cut_2a; s.predecessor[ref_successor_cut_2a] = ref_predecessor_cut_1a;
																	s.predecessor[ref_successor_cut_1b] = ref_predecessor_cut_2b; s.predecessor[ref_successor_cut_2b] = ref_predecessor_cut_1b;
																	s.successor[ref_predecessor_cut_1a] = ref_successor_cut_2a; s.successor[ref_predecessor_cut_2a] = ref_successor_cut_1a;
																	s.successor[ref_predecessor_cut_1b] = ref_successor_cut_2b; s.successor[ref_predecessor_cut_2b] = ref_successor_cut_1b;

																	// Check the feasibility of the routes involved:
																	if (check_load(p, s, ref_predecessor_cut_1a, ref_successor_cut_1b) == true && check_service_start_Tang(p, s, vehicle_1) == true
																		&& check_load(p, s,  ref_predecessor_cut_2a, ref_successor_cut_2b) == true && check_service_start_Tang(p, s, vehicle_2) == true
																		&& (transfers == 0 || check_service_start_BellmanFord(p, s) == true)) {
																		
																		best_exchange_cost = exchange_cost;
																		best_vehicle_1 = vehicle_1; best_vehicle_2 = vehicle_2;
																		best_predecessor_cut_1a = ref_predecessor_cut_1a; best_predecessor_cut_2a = ref_predecessor_cut_2a;
																		best_predecessor_cut_1b = ref_predecessor_cut_1b; best_predecessor_cut_2b = ref_predecessor_cut_2b;
																		best_successor_cut_1a = ref_successor_cut_1a; best_successor_cut_2a = ref_successor_cut_2a;
																		best_successor_cut_1b = ref_successor_cut_1b; best_successor_cut_2b = ref_successor_cut_2b;
																	}

																	// Undo the exchange:
																	s.predecessor[ref_successor_cut_1a] = ref_predecessor_cut_1a; s.predecessor[ref_successor_cut_2a] = ref_predecessor_cut_2a;
																	s.predecessor[ref_successor_cut_1b] = ref_predecessor_cut_1b; s.predecessor[ref_successor_cut_2b] = ref_predecessor_cut_2b;
																	s.successor[ref_predecessor_cut_1a] = ref_successor_cut_1a; s.successor[ref_predecessor_cut_2a] = ref_successor_cut_2a;
																	s.successor[ref_predecessor_cut_1b] = ref_successor_cut_1b; s.successor[ref_predecessor_cut_2b] = ref_successor_cut_2b;

																	ref = ref_successor_cut_1a;
																	while (true) {
																		s.map_to_vehicle[ref] = vehicle_1;
																		ref = s.successor[ref];
																		if (ref == ref_successor_cut_1b) { break; }
																	}

																	ref = ref_successor_cut_2a;
																	while (true) {
																		s.map_to_vehicle[ref] = vehicle_2;
																		ref = s.successor[ref];
																		if (ref == ref_successor_cut_2b) { break; }
																	}
																}
															}
														}
														ref_predecessor_cut_2b = s.successor[ref_predecessor_cut_2b];
														node_predecessor_cut_2b = s.map_to_node[ref_predecessor_cut_2b];
														ref_successor_cut_2b = s.successor[ref_predecessor_cut_2b];
														node_successor_cut_2b = s.map_to_node[ref_successor_cut_2b];
														if (ref_predecessor_cut_2b == ref_end_depot_2) { break; }
													}
												}
											}
											ref_predecessor_cut_2a = s.successor[ref_predecessor_cut_2a];
											node_predecessor_cut_2a = s.map_to_node[ref_predecessor_cut_2a];
											ref_successor_cut_2a = s.successor[ref_predecessor_cut_2a];
											node_successor_cut_2a = s.map_to_node[ref_successor_cut_2a];
											if (ref_predecessor_cut_2a == s.predecessor[ref_end_depot_2]) { break; }
										}
									}
								}
							}
							ref_predecessor_cut_1b = s.successor[ref_predecessor_cut_1b];
							node_predecessor_cut_1b = s.map_to_node[ref_predecessor_cut_1b];
							ref_successor_cut_1b = s.successor[ref_predecessor_cut_1b];
							node_successor_cut_1b = s.map_to_node[ref_successor_cut_1b];
							if (ref_predecessor_cut_1b == ref_end_depot_1) { break; }
						}
					}
					ref_predecessor_cut_1a = s.successor[ref_predecessor_cut_1a];
					node_predecessor_cut_1a = s.map_to_node[ref_predecessor_cut_1a];
					ref_successor_cut_1a = s.successor[ref_predecessor_cut_1a];
					node_successor_cut_1a = s.map_to_node[ref_successor_cut_1a];
					if (ref_predecessor_cut_1a == s.predecessor[ref_end_depot_1]) { break; }
				}
			}
		}

		// Perform the best exchange found (if any):
		if (best_exchange_cost >= -0.001) { break; }
		
		s.total_distance = s.total_distance + best_exchange_cost;

		s.predecessor[best_successor_cut_1a] = best_predecessor_cut_2a; s.predecessor[best_successor_cut_2a] = best_predecessor_cut_1a;
		s.predecessor[best_successor_cut_1b] = best_predecessor_cut_2b; s.predecessor[best_successor_cut_2b] = best_predecessor_cut_1b;
		s.successor[best_predecessor_cut_1a] = best_successor_cut_2a; s.successor[best_predecessor_cut_2a] = best_successor_cut_1a;
		s.successor[best_predecessor_cut_1b] = best_successor_cut_2b; s.successor[best_predecessor_cut_2b] = best_successor_cut_1b;

		// Update all route assignments that have changed (between both cuts):
		int ref = best_successor_cut_1a;
		while (true) {
			s.map_to_vehicle[ref] = best_vehicle_2;
			ref = s.successor[ref];
			if (ref == best_successor_cut_2b) { break; }
		}

		ref = best_successor_cut_2a;
		while (true) {
			s.map_to_vehicle[ref] = best_vehicle_1;
			ref = s.successor[ref];
			if (ref == best_successor_cut_1b) { break; }
		}

		update_load(p, s, best_predecessor_cut_1a, best_successor_cut_1b);
		update_load(p, s, best_predecessor_cut_2a, best_successor_cut_2b);
		update_earliest_time(p, s, best_predecessor_cut_1a);
		update_earliest_time(p, s, best_predecessor_cut_2a);
		update_latest_time(p, s, best_successor_cut_1b);
		update_latest_time(p, s, best_successor_cut_2b);
	}
}


void remove_request(struct problem &p, struct solution &s, int request) {

	int pickup = request;
	int delivery = request + p.n_requests;

	// Case 1 - The request is not partly served by public transport:
	if (transfers == 0 || s.map_to_vehicle[request + 2 * p.n_requests + 2 * p.n_vehicles] == -1) {

		// Case 1a - The pickup and delivery node are served right after each other:
		if (delivery == s.successor[pickup]) {

			int ref_predecessor_pickup = s.predecessor[pickup];
			int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
			int ref_successor_delivery = s.successor[delivery];
			int node_successor_delivery = s.map_to_node[ref_successor_delivery];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_pickup * p.n_nodes + pickup] - p.matrix[pickup * p.n_nodes + delivery]
				- p.matrix[delivery * p.n_nodes + node_successor_delivery] + p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_delivery];

			s.predecessor[ref_successor_delivery] = ref_predecessor_pickup;
			s.successor[ref_predecessor_pickup] = ref_successor_delivery;
		}

		// Case 1b - The pickup and delivery node are not served right after each other:
		else {

			int ref_predecessor_pickup = s.predecessor[pickup];
			int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
			int ref_successor_pickup = s.successor[pickup];
			int node_successor_pickup = s.map_to_node[ref_successor_pickup];
			int ref_predecessor_delivery = s.predecessor[delivery];
			int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
			int ref_successor_delivery = s.successor[delivery];
			int node_successor_delivery = s.map_to_node[ref_successor_delivery];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_pickup * p.n_nodes + pickup] - p.matrix[pickup * p.n_nodes + node_successor_pickup] - p.matrix[node_predecessor_delivery * p.n_nodes + delivery]
				- p.matrix[delivery * p.n_nodes + node_successor_delivery] + p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup] + p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];

			s.predecessor[ref_successor_pickup] = ref_predecessor_pickup;
			s.successor[ref_predecessor_pickup] = ref_successor_pickup;
			s.predecessor[ref_successor_delivery] = ref_predecessor_delivery;
			s.successor[ref_predecessor_delivery] = ref_successor_delivery;
		}

		// Change solution attributes without updating load, earliest/latest time and service start of other nodes (solution remains feasible):
		s.map_to_vehicle[pickup] = -1; s.map_to_vehicle[delivery] = -1;
		s.predecessor[pickup] = -1; s.successor[pickup] = -1;
		s.predecessor[delivery] = -1; s.successor[delivery] = -1;

		s.load_mobile[pickup] = 0; s.load_mobile[delivery] = 0;
		s.load_wheelchair[pickup] = 0; s.load_wheelchair[delivery] = 0;

		s.service_start[pickup] = 0.0; s.service_start[delivery] = 0.0;
		s.earliest_time[pickup] = 0.0; s.latest_time[pickup] = 1440.0;
		s.earliest_time[delivery] = 0.0; s.latest_time[delivery] = 1440.0;
	}

	// Case 2 - The request is partly served by public transport:
	else {

		int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request;
		int node_terminal1 = s.map_to_node[ref_terminal1];
		int ref_terminal2 = 3 * p.n_requests + 2 * p.n_vehicles + request;
		int node_terminal2 = s.map_to_node[ref_terminal2];

		// Case 2a - The pickup and transfer node are served right after each other:
		if (ref_terminal1 == s.successor[pickup]) {

			int ref_predecessor_pickup = s.predecessor[pickup];
			int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
			int ref_successor_terminal1 = s.successor[ref_terminal1];
			int node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_pickup * p.n_nodes + pickup] - p.matrix[pickup * p.n_nodes + node_terminal1]
				- p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1] + p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_terminal1];

			s.predecessor[s.successor[ref_terminal1]] = s.predecessor[pickup];
			s.successor[s.predecessor[pickup]] = s.successor[ref_terminal1];
		}

		// Case 2b - The pickup and transfer node are not served right after each other:
		else {

			int ref_predecessor_pickup = s.predecessor[pickup];
			int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
			int ref_successor_pickup = s.successor[pickup];
			int node_successor_pickup = s.map_to_node[ref_successor_pickup];
			int ref_predecessor_terminal1 = s.predecessor[ref_terminal1];
			int node_predecessor_terminal1 = s.map_to_node[ref_predecessor_terminal1];
			int ref_successor_terminal1 = s.successor[ref_terminal1];
			int node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_pickup * p.n_nodes + pickup] - p.matrix[pickup * p.n_nodes + node_successor_pickup] - p.matrix[node_predecessor_terminal1 * p.n_nodes + node_terminal1]
				- p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1] + p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup] + p.matrix[node_predecessor_terminal1 * p.n_nodes + node_successor_terminal1];

			s.predecessor[s.successor[pickup]] = s.predecessor[pickup];
			s.successor[s.predecessor[pickup]] = s.successor[pickup];
			s.predecessor[s.successor[ref_terminal1]] = s.predecessor[ref_terminal1];
			s.successor[s.predecessor[ref_terminal1]] = s.successor[ref_terminal1];
		}

		// Case 2c - The transfer node and delivery node are served right after each other:
		if (delivery == s.successor[ref_terminal2]) {

			int ref_predecessor_terminal2 = s.predecessor[ref_terminal2];
			int node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];
			int ref_successor_delivery = s.successor[delivery];
			int node_successor_delivery = s.map_to_node[ref_successor_delivery];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] - p.matrix[node_terminal2 * p.n_nodes + delivery]
				- p.matrix[delivery * p.n_nodes + node_successor_delivery] + p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_delivery];

			s.predecessor[s.successor[delivery]] = s.predecessor[ref_terminal2];
			s.successor[s.predecessor[ref_terminal2]] = s.successor[delivery];
		}

		// Case 2d - The transfer node and delivery node are not served right after each other:
		else {

			int ref_predecessor_terminal2 = s.predecessor[ref_terminal2];
			int node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];
			int ref_successor_terminal2 = s.successor[ref_terminal2];
			int node_successor_terminal2 = s.map_to_node[ref_successor_terminal2];
			int ref_predecessor_delivery = s.predecessor[delivery];
			int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
			int ref_successor_delivery = s.successor[delivery];
			int node_successor_delivery = s.map_to_node[ref_successor_delivery];

			s.total_distance = s.total_distance - p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] - p.matrix[node_terminal2 * p.n_nodes + node_successor_terminal2] - p.matrix[node_predecessor_delivery * p.n_nodes + delivery]
				- p.matrix[delivery * p.n_nodes + node_successor_delivery] + p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_terminal2] + p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];

			s.predecessor[s.successor[ref_terminal2]] = s.predecessor[ref_terminal2];
			s.successor[s.predecessor[ref_terminal2]] = s.successor[ref_terminal2];
			s.predecessor[s.successor[delivery]] = s.predecessor[delivery];
			s.successor[s.predecessor[delivery]] = s.successor[delivery];
		}

		// Change solution attributes without updating load, earliest/latest time and service start of other nodes (solution remains feasible):
		s.map_to_vehicle[pickup] = -1; s.map_to_vehicle[delivery] = -1;
		s.map_to_vehicle[ref_terminal1] = -1; s.map_to_vehicle[ref_terminal2] = -1;
		s.map_to_node[ref_terminal1] = -1; s.map_to_node[ref_terminal2] = -1;

		s.predecessor[pickup] = -1; s.successor[pickup] = -1;
		s.predecessor[ref_terminal1] = -1; s.successor[ref_terminal1] = -1;
		s.predecessor[ref_terminal2] = -1; s.successor[ref_terminal2] = -1;
		s.predecessor[delivery] = -1; s.successor[delivery] = -1;

		s.load_mobile[pickup] = 0; s.load_mobile[delivery] = 0;
		s.load_mobile[ref_terminal1] = 0; s.load_mobile[ref_terminal2] = 0;
		s.load_wheelchair[pickup] = 0; s.load_wheelchair[delivery] = 0;
		s.load_wheelchair[ref_terminal1] = 0; s.load_wheelchair[ref_terminal2] = 0;

		s.service_start[pickup] = 0.0; s.service_start[delivery] = 0.0;
		s.service_start[ref_terminal1] = 0.0; s.service_start[ref_terminal2] = 0.0;
		s.earliest_time[pickup] = 0.0; s.latest_time[pickup] = 1440.0;
		s.earliest_time[ref_terminal1] = 0.0; s.latest_time[ref_terminal1] = 1440.0;
		s.earliest_time[ref_terminal2] = 0.0; s.latest_time[ref_terminal2] = 1440.0;
		s.earliest_time[delivery] = 0.0; s.latest_time[delivery] = 1440.0;
	}
}


int clear_terminal(struct problem &p, struct solution &s, int request1_id) {

	int number_removed = 0;

	int ref_request1_terminal1 = request1_id + 2 * p.n_requests + 2 * p.n_vehicles;
	int ref_request1_terminal2 = ref_request1_terminal1 + p.n_requests;

	if (s.map_to_vehicle[ref_request1_terminal1] != -1) {

		for (int request2_id = 0; request2_id < p.n_requests; request2_id++) {

			int ref_request2_terminal1 = request2_id + 2 * p.n_requests + 2 * p.n_vehicles;
			int ref_request2_terminal2 = ref_request2_terminal1 + p.n_requests;

			if (request2_id != request1_id && s.map_to_vehicle[ref_request2_terminal1] != -1 &&
				((s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal1] + 30) ||
				(s.map_to_node[ref_request1_terminal1] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal1] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal1] + 30) ||
				(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal1] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal1] && s.service_start[ref_request2_terminal1] <= s.service_start[ref_request1_terminal2] + 30) ||
				(s.map_to_node[ref_request1_terminal2] == s.map_to_node[ref_request2_terminal2] && s.service_start[ref_request1_terminal2] - 30 <= s.service_start[ref_request2_terminal2] && s.service_start[ref_request2_terminal2] <= s.service_start[ref_request1_terminal2] + 30))) {

				remove_request(p, s, request2_id);
				number_removed += 1;
			}
		}
	}

	return number_removed;
}


void insert_request(struct problem &p, struct solution &s, int request, int vehicle) {

	int pickup = request;
	int delivery = request + p.n_requests;

	s.map_to_vehicle[pickup] = vehicle;
	s.map_to_vehicle[delivery] = vehicle;

	s.total_distance = s.total_distance + s.best_insertion_cost[request];

	s.predecessor[pickup] = s.best_predecessor[pickup];
	s.successor[pickup] = s.best_successor[pickup];
	s.successor[s.predecessor[pickup]] = pickup;
	s.predecessor[s.successor[pickup]] = pickup;

	s.predecessor[delivery] = s.best_predecessor[delivery];
	s.successor[delivery] = s.best_successor[delivery];
	s.successor[s.predecessor[delivery]] = delivery;
	s.predecessor[s.successor[delivery]] = delivery;

	update_load(p, s, pickup, delivery);
	update_earliest_time(p, s, pickup);
	update_latest_time(p, s, delivery);
}


void insert_request_with_transfer(struct problem &p, struct solution &s, int request, int vehicle1, int vehicle2) {

	int pickup = request;
	int delivery = request + p.n_requests;
	int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request;
	int ref_terminal2 = 3 * p.n_requests + 2 * p.n_vehicles + request;

	s.map_to_vehicle[pickup] = vehicle1;
	s.map_to_vehicle[ref_terminal1] = vehicle1;
	s.map_to_vehicle[ref_terminal2] = vehicle2;
	s.map_to_vehicle[delivery] = vehicle2;

	s.map_to_node[ref_terminal1] = s.best_terminal1[request];
	s.map_to_node[ref_terminal2] = s.best_terminal2[request];

	s.total_distance = s.total_distance + s.best_insertion_cost[request];

	s.predecessor[pickup] = s.best_predecessor[pickup];
	s.successor[pickup] = s.best_successor[pickup];
	s.successor[s.predecessor[pickup]] = pickup;
	s.predecessor[s.successor[pickup]] = pickup;

	s.predecessor[ref_terminal1] = s.best_predecessor[ref_terminal1];
	s.successor[ref_terminal1] = s.best_successor[ref_terminal1];
	s.successor[s.predecessor[ref_terminal1]] = ref_terminal1;
	s.predecessor[s.successor[ref_terminal1]] = ref_terminal1;

	s.predecessor[delivery] = s.best_predecessor[delivery];
	s.successor[delivery] = s.best_successor[delivery];
	s.successor[s.predecessor[delivery]] = delivery;
	s.predecessor[s.successor[delivery]] = delivery;

	s.predecessor[ref_terminal2] = s.best_predecessor[ref_terminal2];
	s.successor[ref_terminal2] = s.best_successor[ref_terminal2];
	s.successor[s.predecessor[ref_terminal2]] = ref_terminal2;
	s.predecessor[s.successor[ref_terminal2]] = ref_terminal2;

	update_load(p, s, pickup, ref_terminal1);
	update_load(p, s, ref_terminal2, delivery);

	update_earliest_time(p, s, pickup);
	update_earliest_time(p, s, ref_terminal2);
	update_latest_time(p, s, ref_terminal1);
	update_latest_time(p, s, delivery);
}


void compute_best_insertion_without_pt(struct problem &p, struct solution &s, int request, int vehicle, bool track_second_best) {

	int pickup = request;
	int delivery = request + p.n_requests;

	int ref_begin_depot = 2 * p.n_requests + vehicle;
	int ref_end_depot = 2 * p.n_requests + p.n_vehicles + vehicle;

	double insertion_cost;

	// Loop 1 - Consider all possible positions for the pickup node, starting from the first position after the starting depot:
	int ref_predecessor_pickup = ref_begin_depot;
	int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
	int ref_successor_pickup = s.successor[ref_begin_depot];
	int node_successor_pickup = s.map_to_node[ref_successor_pickup];

	do {
		
		// Preliminary capacity and time window check at the pickup node:
		if (s.earliest_time[ref_predecessor_pickup] + p.nodes[node_predecessor_pickup].service_dur + p.matrix[node_predecessor_pickup * p.n_nodes + pickup] <= p.nodes[pickup].upper_tw
			&& p.nodes[pickup].lower_tw + p.nodes[pickup].service_dur + p.matrix_nontightened[pickup * p.n_nodes + node_successor_pickup] <= s.latest_time[ref_successor_pickup]
			&& s.load_mobile[ref_predecessor_pickup] + p.requests[request].users_mobile <= p.vehicles[vehicle].cap_mobile
			&& s.load_wheelchair[ref_predecessor_pickup] + p.requests[request].users_wheelchair <= p.vehicles[vehicle].cap_wheelchair) {
			
			// Loop 2 - Consider all possible positions for the delivery node, starting from the first position after the pickup node:
			int ref_predecessor_delivery = pickup;
			int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
			int ref_successor_delivery = ref_successor_pickup;
			int node_successor_delivery = s.map_to_node[ref_successor_delivery];

			do {
				
				// Preliminary time window checks at the delivery location:
				if (p.nodes[delivery].lower_tw + p.nodes[delivery].service_dur + p.matrix[delivery * p.n_nodes + node_successor_delivery] <= s.latest_time[ref_successor_delivery]
					&& s.earliest_time[ref_predecessor_delivery] + p.nodes[node_predecessor_delivery].service_dur + p.matrix[node_predecessor_delivery * p.n_nodes + delivery] <= p.nodes[delivery].upper_tw) {
					
					// Compute the insertion cost of the new request:
					if (ref_predecessor_delivery == pickup) {
						insertion_cost = p.matrix[node_predecessor_pickup * p.n_nodes + pickup] + p.matrix[pickup * p.n_nodes + delivery] + p.matrix[delivery * p.n_nodes + node_successor_delivery]
							- p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_delivery];
					}
					else {
						insertion_cost = p.matrix[node_predecessor_pickup * p.n_nodes + pickup] + p.matrix[pickup * p.n_nodes + node_successor_pickup] + p.matrix[node_predecessor_delivery * p.n_nodes + delivery] + p.matrix[delivery * p.n_nodes + node_successor_delivery]
							- p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup] - p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];
					}
					
					// If this is the best-found insertion cost so far, check all feasibility constraints:
					if ((track_second_best == true && insertion_cost < s.second_best_insertion_cost[request] - 0.001) ||
						(track_second_best == false && insertion_cost < s.best_insertion_cost[request] - 0.001)) {
						
						// Temporaly replace the predecessors and successors in the solution:
						s.map_to_vehicle[pickup] = vehicle; s.map_to_vehicle[delivery] = vehicle;
						s.predecessor[pickup] = ref_predecessor_pickup; s.successor[ref_predecessor_pickup] = pickup;
						s.predecessor[ref_successor_pickup] = pickup; s.successor[pickup] = ref_successor_pickup;
						s.predecessor[delivery] = ref_predecessor_delivery; s.successor[ref_predecessor_delivery] = delivery;
						s.predecessor[ref_successor_delivery] = delivery; s.successor[delivery] = ref_successor_delivery;

						// Check feasibility with respect to time and load:
						if (check_load(p, s, pickup, delivery) == true &&
							check_service_start_Tang(p, s, vehicle) == true &&
							(transfers == 0 || check_service_start_BellmanFord(p, s) == true)) {

							if (insertion_cost < s.best_insertion_cost[request] - 0.001) {

								s.second_best_insertion_cost[request] = s.best_insertion_cost[request];
								s.best_insertion_cost[request] = insertion_cost;

								s.second_best_vehicle1[request] = s.best_vehicle1[request];
								s.second_best_vehicle2[request] = s.best_vehicle2[request];
								s.best_vehicle1[request] = vehicle;
								s.best_vehicle2[request] = -1;

								s.best_predecessor[pickup] = ref_predecessor_pickup;
								s.best_successor[pickup] = ref_successor_pickup;
								s.best_predecessor[delivery] = ref_predecessor_delivery;
								s.best_successor[delivery] = ref_successor_delivery;
							}

							else if (track_second_best == true) {
							
								s.second_best_insertion_cost[request] = insertion_cost;

								s.second_best_vehicle1[request] = vehicle;
								s.second_best_vehicle2[request] = -1;
							}
						}

						// Restore the original predecessors and successors in the solution:
						s.map_to_vehicle[pickup] = -1; s.map_to_vehicle[delivery] = -1;
						s.successor[ref_predecessor_delivery] = ref_successor_delivery;
						s.predecessor[ref_successor_delivery] = ref_predecessor_delivery;
						s.successor[ref_predecessor_pickup] = ref_successor_pickup;
						s.predecessor[ref_successor_pickup] = ref_predecessor_pickup;
					}
				}

				// Consider the next position for the delivery node:
				ref_predecessor_delivery = ref_successor_delivery;
				node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];
				ref_successor_delivery = s.successor[ref_predecessor_delivery];
				node_successor_delivery = s.map_to_node[ref_successor_delivery];

				// Break the loop as soon as the time window of the delivery node can no longer be respected:
				if (s.earliest_time[ref_predecessor_delivery] + p.nodes[node_predecessor_delivery].service_dur + p.matrix_nontightened[node_predecessor_delivery * p.n_nodes + delivery] > p.nodes[delivery].upper_tw) { break; }

			} while (ref_successor_delivery != ref_begin_depot);
		}

		// Consider the next position for the pickup node:
		ref_predecessor_pickup = ref_successor_pickup;
		node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
		ref_successor_pickup = s.successor[ref_predecessor_pickup];
		node_successor_pickup = s.map_to_node[ref_successor_pickup];

		// Break the loop as soon as the time window of the pickup node can no longer be respected:
		if (s.earliest_time[ref_predecessor_pickup] + p.nodes[node_predecessor_pickup].service_dur + p.matrix_nontightened[node_predecessor_pickup * p.n_nodes + pickup] > p.nodes[pickup].upper_tw) { break; }

	} while (ref_successor_pickup != ref_begin_depot);
}


void compute_best_insertion_pt(struct problem &p, struct solution &s, int request, int vehicle1, int vehicle2, int terminal1, int terminal2, bool track_second_best) {
	
	int pickup = request;
	int delivery = request + p.n_requests;
	int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request;
	int ref_terminal2 = 3 * p.n_requests + 2 * p.n_vehicles + request;

	int node_terminal1 = 2 * p.n_requests + p.n_vehicles + terminal1;
	int node_terminal2 = 2 * p.n_requests + p.n_vehicles + terminal2;

	int ref_begin_depot1 = 2 * p.n_requests + vehicle1;
	int ref_begin_depot2 = 2 * p.n_requests + vehicle2;
	int ref_end_depot1 = 2 * p.n_requests + p.n_vehicles + vehicle1;
	int ref_end_depot2 = 2 * p.n_requests + p.n_vehicles + vehicle2;

	double insertion_cost;

	// Check whether a trip via the chosen transfers terminals can respect the maximum user ride time:
	if (p.matrix[pickup * p.n_nodes + node_terminal1] + p.nodes[node_terminal1].service_dur + p.matrix_pt[terminal1 * p.n_terminals + terminal2] + p.nodes[node_terminal2].service_dur + p.matrix[node_terminal2 * p.n_nodes + delivery] < p.requests[request].max_ride_time) {

		double latest_arrival_terminal1 = p.nodes[delivery].upper_tw - p.matrix[node_terminal2 * p.n_nodes + delivery] - p.nodes[node_terminal2].service_dur - p.matrix_pt[terminal1 * p.n_terminals + terminal2] - p.nodes[node_terminal1].service_dur;
		double earliest_start_terminal2 = p.nodes[pickup].lower_tw + p.nodes[pickup].service_dur + p.matrix[pickup * p.n_nodes + node_terminal1] + p.nodes[node_terminal1].service_dur + p.matrix_pt[terminal1 * p.n_terminals + terminal2];

		// Loop 1 - Consider all possible positions for the pickup node, starting from the first position after the starting depot:
		int ref_predecessor_pickup = ref_begin_depot1;
		int node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
		int ref_successor_pickup = s.successor[ref_begin_depot1];
		int node_successor_pickup = s.map_to_node[ref_successor_pickup];

		do {

			// Preliminary capacity and time window check at the pickup node:
			if (s.earliest_time[ref_predecessor_pickup] + p.nodes[node_predecessor_pickup].service_dur + p.matrix[node_predecessor_pickup * p.n_nodes + pickup] <= p.nodes[pickup].upper_tw
				&& p.nodes[pickup].lower_tw + p.nodes[pickup].service_dur + p.matrix_nontightened[pickup * p.n_nodes + node_successor_pickup] <= s.latest_time[ref_successor_pickup]	 // Use matrix_nontightened to allow correct arc elimination
				&& s.load_mobile[ref_predecessor_pickup] + p.requests[request].users_mobile <= p.vehicles[vehicle1].cap_mobile
				&& s.load_wheelchair[ref_predecessor_pickup] + p.requests[request].users_wheelchair <= p.vehicles[vehicle1].cap_wheelchair) {

				// Loop 2 - Consider all possible positions for the delivery node in reverse order, starting from the last position before the ending depot:
				int ref_successor_delivery = ref_end_depot2;
				int node_successor_delivery = s.map_to_node[ref_successor_delivery];
				int ref_predecessor_delivery = s.predecessor[ref_end_depot2];
				int node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];

				do {

					// Preliminary time window checks at the delivery location:
					if (p.nodes[delivery].lower_tw + p.nodes[delivery].service_dur + p.matrix[delivery * p.n_nodes + node_successor_delivery] <= s.latest_time[ref_successor_delivery]
						&& s.earliest_time[ref_predecessor_delivery] + p.nodes[node_predecessor_delivery].service_dur + p.matrix_nontightened[node_predecessor_delivery * p.n_nodes + delivery] <= p.nodes[delivery].upper_tw) {	 // Use matrix_nontightened to allow correct arc elimination

						// Loop 3 - Consider all possible positions for the first transfer node, starting from the first position after the pickup node:
						int ref_predecessor_terminal1 = pickup;
						int node_predecessor_terminal1 = s.map_to_node[ref_predecessor_terminal1];
						int ref_successor_terminal1 = ref_successor_pickup;
						int node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

						do {

							// Loop 4 - Consider all possible positions for the second transfer node in reverse order, starting from the first position before the delivery node:
							int ref_successor_terminal2 = delivery;
							int node_successor_terminal2 = s.map_to_node[ref_successor_terminal2];
							int ref_predecessor_terminal2 = ref_predecessor_delivery;
							int node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];

							do {

								// Compute the insertion cost of the new request:
								if (ref_predecessor_terminal1 == pickup) {
									insertion_cost = p.matrix[node_predecessor_pickup * p.n_nodes + pickup] + p.matrix[pickup * p.n_nodes + node_terminal1] + p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1]
										- p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_terminal1];
								}
								else {
									insertion_cost = p.matrix[node_predecessor_pickup * p.n_nodes + pickup] + p.matrix[pickup * p.n_nodes + node_successor_pickup] + p.matrix[node_predecessor_terminal1 * p.n_nodes + node_terminal1] + p.matrix[node_terminal1 * p.n_nodes + node_successor_terminal1]
										- p.matrix[node_predecessor_pickup * p.n_nodes + node_successor_pickup] - p.matrix[node_predecessor_terminal1 * p.n_nodes + node_successor_terminal1];
								}
								if (delivery == ref_successor_terminal2) {
									insertion_cost = insertion_cost + p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] + p.matrix[node_terminal2 * p.n_nodes + delivery] + p.matrix[delivery * p.n_nodes + node_successor_delivery]
										- p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_delivery];
								}
								else {
									insertion_cost = insertion_cost + p.matrix[node_predecessor_terminal2 * p.n_nodes + node_terminal2] + p.matrix[node_terminal2 * p.n_nodes + node_successor_terminal2] + p.matrix[node_predecessor_delivery * p.n_nodes + delivery] + p.matrix[delivery * p.n_nodes + node_successor_delivery]
										- p.matrix[node_predecessor_terminal2 * p.n_nodes + node_successor_terminal2] - p.matrix[node_predecessor_delivery * p.n_nodes + node_successor_delivery];
								}
									
								// If this is the best-found insertion cost so far, check all feasibility constraints:
								if ((track_second_best == true && insertion_cost < s.second_best_insertion_cost[request] - 0.001) ||
									(track_second_best == false && insertion_cost < s.best_insertion_cost[request] - 0.001)) {

									// Temporaly replace the predecessors and successors in the solution:
									s.map_to_vehicle[pickup] = vehicle1;
									s.map_to_vehicle[delivery] = vehicle2;
									s.map_to_vehicle[ref_terminal1] = vehicle1;
									s.map_to_vehicle[ref_terminal2] = vehicle2;
									s.map_to_node[ref_terminal1] = 2 * p.n_requests + p.n_vehicles + terminal1;
									s.map_to_node[ref_terminal2] = 2 * p.n_requests + p.n_vehicles + terminal2;

									if (ref_predecessor_terminal1 == pickup) {
										s.predecessor[pickup] = ref_predecessor_pickup; s.successor[ref_predecessor_pickup] = pickup;
										s.predecessor[ref_terminal1] = pickup; s.successor[pickup] = ref_terminal1;
										s.predecessor[ref_successor_terminal1] = ref_terminal1; s.successor[ref_terminal1] = ref_successor_terminal1;
									}
									else {
										s.predecessor[pickup] = ref_predecessor_pickup; s.successor[ref_predecessor_pickup] = pickup;
										s.predecessor[ref_successor_pickup] = pickup; s.successor[pickup] = ref_successor_pickup;
										s.predecessor[ref_terminal1] = ref_predecessor_terminal1; s.successor[ref_predecessor_terminal1] = ref_terminal1;
										s.predecessor[ref_successor_terminal1] = ref_terminal1; s.successor[ref_terminal1] = ref_successor_terminal1;
									}

									if (delivery == ref_successor_terminal2) {
										s.predecessor[ref_terminal2] = ref_predecessor_terminal2; s.successor[ref_predecessor_terminal2] = ref_terminal2;
										s.predecessor[delivery] = ref_terminal2; s.successor[ref_terminal2] = delivery;
										s.predecessor[ref_successor_delivery] = delivery; s.successor[delivery] = ref_successor_delivery;
									}
									else {
										s.predecessor[ref_terminal2] = ref_predecessor_terminal2; s.successor[ref_predecessor_terminal2] = ref_terminal2;
										s.predecessor[ref_successor_terminal2] = ref_terminal2; s.successor[ref_terminal2] = ref_successor_terminal2;
										s.predecessor[delivery] = ref_predecessor_delivery; s.successor[ref_predecessor_delivery] = delivery;
										s.predecessor[ref_successor_delivery] = delivery; s.successor[delivery] = ref_successor_delivery;
									}

									// Check feasibility with respect to time and load:
									if (check_load(p, s, pickup, ref_terminal1) == true && check_load(p, s, ref_terminal2, delivery) == true &&
										check_service_start_Tang(p, s, vehicle1) == true && check_service_start_Tang(p, s, vehicle2) == true &&
										(transfers == 0 || check_service_start_BellmanFord(p, s) == true)) {

										if (insertion_cost < s.best_insertion_cost[request] - 0.001) {

											s.second_best_insertion_cost[request] = s.best_insertion_cost[request];
											s.best_insertion_cost[request] = insertion_cost;

											s.second_best_vehicle1[request] = s.best_vehicle1[request];
											s.second_best_vehicle2[request] = s.best_vehicle2[request];
											s.best_vehicle1[request] = vehicle1;
											s.best_vehicle2[request] = vehicle2;

											s.best_terminal1[request] = node_terminal1;
											s.best_terminal2[request] = node_terminal2;

											if (ref_predecessor_terminal1 == pickup) {
												s.best_predecessor[pickup] = ref_predecessor_pickup;
												s.best_successor[pickup] = ref_terminal1;
												s.best_predecessor[ref_terminal1] = pickup;
												s.best_successor[ref_terminal1] = ref_successor_terminal1;
											}
											else {
												s.best_predecessor[pickup] = ref_predecessor_pickup;
												s.best_successor[pickup] = ref_successor_pickup;
												s.best_predecessor[ref_terminal1] = ref_predecessor_terminal1;
												s.best_successor[ref_terminal1] = ref_successor_terminal1;
											}

											if (delivery == ref_successor_terminal2) {
												s.best_predecessor[ref_terminal2] = ref_predecessor_terminal2;
												s.best_successor[ref_terminal2] = delivery;
												s.best_predecessor[delivery] = ref_terminal2;
												s.best_successor[delivery] = ref_successor_delivery;
											}
											else {
												s.best_predecessor[ref_terminal2] = ref_predecessor_terminal2;
												s.best_successor[ref_terminal2] = ref_successor_terminal2;
												s.best_predecessor[delivery] = ref_predecessor_delivery;
												s.best_successor[delivery] = ref_successor_delivery;
											}
										}

										else if (track_second_best == true) {
										
											s.second_best_insertion_cost[request] = insertion_cost;

											s.second_best_vehicle1[request] = vehicle1;
											s.second_best_vehicle2[request] = vehicle2;
										}
									}

									// Restore the original predecessors and successors in the solution:
									s.map_to_vehicle[pickup] = -1;
									s.map_to_vehicle[delivery] = -1;
									s.map_to_vehicle[ref_terminal1] = -1;
									s.map_to_vehicle[ref_terminal2] = -1;
									s.map_to_node[ref_terminal1] = -1;
									s.map_to_node[ref_terminal2] = -1;

									if (ref_predecessor_terminal1 == pickup) {
										s.predecessor[pickup] = -1; s.successor[ref_predecessor_pickup] = ref_successor_terminal1;
										s.successor[pickup] = -1; s.predecessor[ref_terminal1] = -1;
										s.predecessor[ref_successor_terminal1] = ref_predecessor_pickup; s.successor[ref_terminal1] = -1;
									}
									else {
										s.predecessor[pickup] = -1; s.successor[ref_predecessor_pickup] = ref_successor_pickup;
										s.predecessor[ref_successor_pickup] = ref_predecessor_pickup; s.successor[pickup] = -1;
										s.predecessor[ref_terminal1] = -1; s.successor[ref_predecessor_terminal1] = ref_successor_terminal1;
										s.predecessor[ref_successor_terminal1] = ref_predecessor_terminal1; s.successor[ref_terminal1] = -1;
									}

									if (delivery == ref_successor_terminal2) {
										s.predecessor[ref_terminal2] = -1; s.successor[ref_predecessor_terminal2] = ref_successor_delivery;
										s.successor[ref_terminal2] = -1; s.predecessor[delivery] = -1;
										s.predecessor[ref_successor_delivery] = ref_predecessor_terminal2; s.successor[delivery] = -1;
									}
									else {
										s.predecessor[ref_terminal2] = -1; s.successor[ref_predecessor_terminal2] = ref_successor_terminal2;
										s.predecessor[ref_successor_terminal2] = ref_predecessor_terminal2; s.successor[ref_terminal2] = -1;
										s.predecessor[delivery] = -1; s.successor[ref_predecessor_delivery] = ref_successor_delivery;
										s.predecessor[ref_successor_delivery] = ref_predecessor_delivery; s.successor[delivery] = -1;
									}
								}

								// Consider the previous position for the second transfer node:
								ref_successor_terminal2 = ref_predecessor_terminal2;
								node_successor_terminal2 = s.map_to_node[ref_successor_terminal2];
								ref_predecessor_terminal2 = s.predecessor[ref_successor_terminal2];
								node_predecessor_terminal2 = s.map_to_node[ref_predecessor_terminal2];

								// Break the loop as soon as the second transfer terminal should be left before the earliest possible time:
								if (s.latest_time[ref_successor_terminal2] - p.matrix_nontightened[node_terminal2 * p.n_nodes + node_successor_terminal2] - p.nodes[node_terminal2].service_dur < earliest_start_terminal2) {
									break;
								}

							} while (ref_predecessor_terminal2 != ref_end_depot2);

							// Consider the next position for the first transfer node:
							ref_predecessor_terminal1 = ref_successor_terminal1;
							node_predecessor_terminal1 = s.map_to_node[ref_predecessor_terminal1];
							ref_successor_terminal1 = s.successor[ref_predecessor_terminal1];
							node_successor_terminal1 = s.map_to_node[ref_successor_terminal1];

							// Break the loop as soon as the first transfer terminal cannot be reached in time:
							if (s.earliest_time[ref_predecessor_terminal1] + p.nodes[node_predecessor_terminal1].service_dur + p.matrix_nontightened[node_predecessor_terminal1 * p.n_nodes + node_terminal1] > latest_arrival_terminal1) {
								break;
							}

						} while (ref_successor_terminal1 != ref_begin_depot1);
					}

					// Consider the previous position for the delivery node:
					ref_successor_delivery = ref_predecessor_delivery;
					node_successor_delivery = s.map_to_node[ref_successor_delivery];
					ref_predecessor_delivery = s.predecessor[ref_successor_delivery];
					node_predecessor_delivery = s.map_to_node[ref_predecessor_delivery];

					// Break the loop as soon as the time window of the delivery node can no longer be respected or the second transfer terminal cannot be left in time:
					if (s.latest_time[ref_successor_delivery] - p.matrix_nontightened[delivery * p.n_nodes + node_successor_delivery] - p.nodes[delivery].service_dur < p.nodes[delivery].lower_tw ||
						s.latest_time[ref_successor_delivery] - p.matrix_nontightened[delivery * p.n_nodes + node_successor_delivery] - p.nodes[delivery].service_dur - p.matrix_nontightened[node_terminal2 * p.n_nodes + delivery] - p.nodes[node_terminal2].service_dur < earliest_start_terminal2) {
						break;
					}

				} while (ref_predecessor_delivery != ref_end_depot2);
			}

			// Consider the next position for the pickup node:
			ref_predecessor_pickup = ref_successor_pickup;
			node_predecessor_pickup = s.map_to_node[ref_predecessor_pickup];
			ref_successor_pickup = s.successor[ref_predecessor_pickup];
			node_successor_pickup = s.map_to_node[ref_successor_pickup];

			// Break the loop as soon as the time window of the pickup node can no longer be respected or the first transfer terminal cannot be reached in time:
			if (s.earliest_time[ref_predecessor_pickup] + p.nodes[node_predecessor_pickup].service_dur + p.matrix_nontightened[node_predecessor_pickup * p.n_nodes + pickup] > p.nodes[pickup].upper_tw ||
				s.earliest_time[ref_predecessor_pickup] + p.nodes[node_predecessor_pickup].service_dur + p.matrix_nontightened[node_predecessor_pickup * p.n_nodes + pickup] + p.nodes[pickup].service_dur + p.matrix_nontightened[pickup * p.n_nodes + node_terminal1] > latest_arrival_terminal1) {
				break;
			}

		} while (ref_successor_pickup != ref_begin_depot1);
	}
}


bool check_load(struct problem &p, struct solution &s, int ref_first_node, int ref_last_node) {

	int ref = ref_first_node;
	int previous_ref = s.predecessor[ref];

	int vehicle = s.map_to_vehicle[ref];
	int start_depot = 2 * p.n_requests + vehicle;
	int end_depot = start_depot + p.n_vehicles;

	// No check needed for the start depot (if included):
	if (ref == start_depot) { previous_ref = ref; ref = s.successor[ref]; }

	int load_mobile = s.load_mobile[previous_ref];
	int load_wheelchair = s.load_wheelchair[previous_ref];

	while (ref != s.successor[ref_last_node]) {

		int request = s.map_to_request[ref];

		// Nodes where users are picked up:
		if (ref < p.n_requests || ref >= 3 * p.n_requests + 2 * p.n_vehicles) {
			load_mobile = load_mobile + p.requests[request].users_mobile;
			load_wheelchair = load_wheelchair + p.requests[request].users_wheelchair;

			// Check vehicle capacity (without upgrading conditions):
			if ((load_mobile > p.vehicles[vehicle].cap_mobile) ||
				(load_wheelchair > p.vehicles[vehicle].cap_wheelchair)) {
				return false;
			}
		}

		// Nodes where users are delivered:
		else if (ref < 2 * p.n_requests || ref >= 2 * p.n_requests + 2 * p.n_vehicles) {
			load_mobile = load_mobile - p.requests[request].users_mobile;
			load_wheelchair = load_wheelchair - p.requests[request].users_wheelchair;
		}

		ref = s.successor[ref];
	}

	return true;
}


bool check_completeness(struct problem &p, struct solution &s) {

	for (int request = 0; request < p.n_requests; request++) {
		if (s.map_to_vehicle[request] == -1 || s.map_to_vehicle[request + p.n_requests] == -1) { return false; }
	}

	return true;
}


bool check_service_start_Tang(struct problem &p, struct solution &s, int vehicle, bool update) {

	double *arrival_time = new double[p.n_references];
	double *service_start = new double[p.n_references];

	int ref_begin_depot = 2 * p.n_requests + vehicle;
	int ref_end_depot = 2 * p.n_requests + p.n_vehicles + vehicle;
	
	// Pass 1 (forward) - Set earliest arrival time and service start without taking into account ride times:

	service_start[ref_begin_depot] = p.nodes[ref_begin_depot].lower_tw;

	int ref_previous_node = ref_begin_depot;
	int node_previous_node = s.map_to_node[ref_begin_depot];
	int ref_current_node = s.successor[ref_begin_depot];
	int node_current_node = s.map_to_node[ref_current_node];

	while (true) {

		arrival_time[ref_current_node] = service_start[ref_previous_node] + p.nodes[node_previous_node].service_dur + p.matrix[node_previous_node * p.n_nodes + node_current_node];
		if (arrival_time[ref_current_node] > p.nodes[node_current_node].upper_tw + 0.001) { delete[] arrival_time; delete[] service_start; return false; }

		service_start[ref_current_node] = max(arrival_time[ref_current_node], p.nodes[node_current_node].lower_tw);

		if (ref_current_node == ref_end_depot) { break; }

		ref_previous_node = ref_current_node;
		node_previous_node = s.map_to_node[ref_previous_node];
		ref_current_node = s.successor[ref_current_node];
		node_current_node = s.map_to_node[ref_current_node];
	}

	// Pass 2 (backward) - Correct ride time violations by postponing the start of service in pickup nodes and adapt the start of service in all following nodes:

	double cumulative_waiting = service_start[ref_current_node] - arrival_time[ref_current_node];

	ref_current_node = s.predecessor[ref_current_node];
	node_current_node = s.map_to_node[ref_current_node];

	while (true) {

		if (p.nodes[node_current_node].type == 0 && s.map_to_vehicle[ref_current_node + p.n_requests] == vehicle) {

			// Compute potential ride time violations:
			double delta = service_start[ref_current_node + p.n_requests] - (service_start[ref_current_node] + p.nodes[node_current_node].service_dur) - p.requests[ref_current_node].max_ride_time;

			if (delta > 0.001) {

				// Terminate the procedure if the ride time violation is irreparable (note: use of "cumulative waiting" speeds up the procedure, but is not necessary):
				if (delta > cumulative_waiting + 0.001) { delete[] arrival_time; delete[] service_start; return false; }

				service_start[ref_current_node] = service_start[ref_current_node] + delta;
				if (service_start[ref_current_node] > p.nodes[node_current_node].upper_tw + 0.001) { delete[] arrival_time; delete[] service_start; return false; }
				cumulative_waiting = cumulative_waiting - delta;

				// Update the start of service in all following nodes:
				int ref_current_node_2 = s.successor[ref_current_node];
				int node_current_node_2 = s.map_to_node[ref_current_node_2];
				int ref_previous_node_2 = ref_current_node;
				int node_previous_node_2 = s.map_to_node[ref_previous_node_2];

				while (true) {

					service_start[ref_current_node_2] = max(service_start[ref_current_node_2], service_start[ref_previous_node_2] + p.nodes[node_previous_node_2].service_dur + p.matrix[node_previous_node_2 * p.n_nodes + node_current_node_2]);
					if (ref_current_node_2 == ref_end_depot) { break; }

					ref_previous_node_2 = ref_current_node_2;
					node_previous_node_2 = s.map_to_node[ref_previous_node_2];
					ref_current_node_2 = s.successor[ref_current_node_2];
					node_current_node_2 = s.map_to_node[ref_current_node_2];
				}
			}
		}

		cumulative_waiting = cumulative_waiting + service_start[ref_current_node] - arrival_time[ref_current_node];

		// Check the maximum route duration constraint (not included in original paper):
		if (ref_current_node == ref_begin_depot) {

			double route_duration = service_start[ref_end_depot] - service_start[ref_begin_depot];
			if (route_duration > p.vehicles[vehicle].max_route_dur + 0.001) { service_start[ref_begin_depot] += (route_duration - p.vehicles[vehicle].max_route_dur); }
			break;
		}

		ref_current_node = s.predecessor[ref_current_node];
		node_current_node = s.map_to_node[ref_current_node];
	}

	// Pass 3 (forward) - Check whether the resulting schedule is feasible:

	ref_previous_node = ref_begin_depot;
	node_previous_node = s.map_to_node[ref_previous_node];
	ref_current_node = s.successor[ref_previous_node];
	node_current_node = s.map_to_node[ref_current_node];

	while (true) {

		arrival_time[ref_current_node] = service_start[ref_previous_node] + p.nodes[node_previous_node].service_dur + p.matrix[node_previous_node * p.n_nodes + node_current_node];
		service_start[ref_current_node] = max(service_start[ref_current_node], arrival_time[ref_current_node]);

		if (service_start[ref_current_node] > p.nodes[node_current_node].upper_tw + 0.001) { delete[] arrival_time; delete[] service_start; return false; }

		if (p.nodes[node_current_node].type == 1 && s.map_to_vehicle[ref_current_node - p.n_requests] == vehicle) {

			// Compute potential ride time violations:
			double delta = service_start[ref_current_node] - (service_start[ref_current_node - p.n_requests] + p.nodes[node_current_node - p.n_requests].service_dur) - p.requests[ref_current_node - p.n_requests].max_ride_time;
			if (delta > 0.001) { delete[] arrival_time; delete[] service_start; return false; }
		}

		if (ref_current_node == ref_end_depot) { break; }
		
		ref_previous_node = ref_current_node;
		node_previous_node = s.map_to_node[ref_previous_node];
		ref_current_node = s.successor[ref_current_node];
		node_current_node = s.map_to_node[ref_current_node];
	}

	// Update the time schedule whenever applicable:
	if (update == true) {
		for (int ref = 0; ref < p.n_references; ref++) {
			if (s.map_to_vehicle[ref] == vehicle) {
				s.service_start[ref] = service_start[ref];
			}
		}
	}

	delete[] arrival_time;
	delete[] service_start;

	return true;
}


bool check_service_start_BellmanFord(struct problem &p, struct solution &s, bool update) {
	
	int *sources = new int[p.n_references + 2 * p.n_requests];
	int *destinations = new int[p.n_references + 2 * p.n_requests];
	double *weights = new double[p.n_references + 2 * p.n_requests];
	double *schedule = new double[p.n_references];

	bool change;
	int count_edges = 0;
	int count_vertices = 0;
	
	// Loop over all (real) references in the network:
	for (int i = 0; i < p.n_references; i++) {

		// Add edges between nodes that are included in the solution:
		if (s.map_to_vehicle[i] != -1) {

			int node = s.map_to_node[i];

			// Set the initial time of service equal to the upper bound of the time window:
			schedule[i] = p.nodes[node].upper_tw;

			count_vertices = count_vertices + 1;

			// Travel times between successive nodes (not for end depots):
			if (i < 2 * p.n_requests + p.n_vehicles || i >= 2 * p.n_requests + 2 * p.n_vehicles) {
				sources[count_edges] = s.successor[i];
				destinations[count_edges] = i;
				weights[count_edges] = -p.nodes[node].service_dur - p.matrix[node * p.n_nodes + s.map_to_node[s.successor[i]]];
				count_edges = count_edges + 1;
			}

			// Maximum route duration (only for start depots):
			if (i >= 2 * p.n_requests && i < 2 * p.n_requests + p.n_vehicles) {
				sources[count_edges] = i;
				destinations[count_edges] = i + p.n_vehicles;
				weights[count_edges] = p.vehicles[s.map_to_vehicle[i]].max_route_dur;
				count_edges = count_edges + 1;
			}

			// Maximum user ride time (only for pickup nodes):
			if (i < p.n_requests) {
				sources[count_edges] = i;
				destinations[count_edges] = i + p.n_requests;
				weights[count_edges] = p.nodes[node].service_dur + p.requests[i].max_ride_time;
				count_edges = count_edges + 1;
			}

			// Synchronisation between routes (only for first transfer nodes):
			if (i >= 2 * p.n_requests + 2 * p.n_vehicles && i < 3 * p.n_requests + 2 * p.n_vehicles) {
				sources[count_edges] = i + p.n_requests;
				destinations[count_edges] = i;
				weights[count_edges] = 0.0;
				count_edges = count_edges + 1;
			}
		}
	}

	// Iteratively search for the lowest feasible time of service (V iterations);
	for (int i = 0; i < count_vertices; i++) {

		change = false;

		// For all arcs, check whether the time of service in the delivery node can be advanced:
		for (int j = 0; j < count_edges; j++) {

			int u = sources[j];
			int v = destinations[j];
			double weight = weights[j];

			// Situation 1 - arcs between two transfer nodes:
			if (u >= 3 * p.n_requests + 2 * p.n_vehicles && v == u - p.n_requests) {

				// Given the time of service at the second transfer node (u), compute how much time has passed since the departure by PT at the first transfer node (v): 
				double required_time_pt = compute_required_time_pt(p, s.map_to_node[v], s.map_to_node[u], schedule[u]);

				if (schedule[u] - required_time_pt < schedule[v] - 0.0001) {
					schedule[v] = schedule[u] - required_time_pt;

					change = true;

					// Declare infeasibility if there is a negative cycle in the solution:
					if (schedule[v] < p.nodes[s.map_to_node[v]].lower_tw - 0.0001) { delete[] sources; delete[] destinations; delete[] weights; delete[] schedule; return false; }
				}
			}

			// Situation 2 - other arcs:
			else if (schedule[u] + weight < schedule[v] - 0.0001) {
				schedule[v] = schedule[u] + weight;

				change = true;

				// Declare infeasibility if there is a negative cycle in the solution:
				if (schedule[v] < p.nodes[s.map_to_node[v]].lower_tw - 0.0001) { delete[] sources; delete[] destinations; delete[] weights; delete[] schedule; return false; }
			}

		}

		// Leave this loop as soon as the solution is stable:
		if (change == false) { break; }
	}

	// Declare infeasibility if the solution is not stable after V iterations:
	if (change == true) { delete[] sources; delete[] destinations; delete[] weights; delete[] schedule; return false; }

	// Update the time schedule whenever applicable:
	if (update == true) {
		for (int ref = 0; ref < p.n_references; ref++) {
			if (s.map_to_vehicle[ref] != -1) {
				s.service_start[ref] = schedule[ref];
			}
		}
	}

	delete[] sources;
	delete[] destinations;
	delete[] weights;
	delete[] schedule;
	return true;
}


double compute_required_time_pt(struct problem &p, int node_origin_terminal, int node_destination_terminal, double pickup_time_second_transfer) {

	double travel_time, service_duration, waiting_time;

	int origin_terminal = node_origin_terminal - 2 * p.n_requests - p.n_vehicles;
	int destination_terminal = node_destination_terminal - 2 * p.n_requests - p.n_vehicles;

	////////////////////////////////////////////////////////////////////
	// CODE BELOW FOR ARTIFICIAL INSTANCES OF MOLENBRUCH ET AL (2018) //
	////////////////////////////////////////////////////////////////////

	// Compute the required travel time:
	travel_time = p.matrix_pt[origin_terminal * p.n_terminals + destination_terminal];

	// Add the service duration at the first transfer terminal:
	service_duration = p.nodes[node_origin_terminal].service_dur;

	// Add the waiting time at the second transfer terminal (depending on which service is taken towards that terminal):
	if (p.nodes[node_destination_terminal].x_coord < p.nodes[node_origin_terminal].x_coord ||
		p.nodes[node_destination_terminal].y_coord > p.nodes[node_origin_terminal].y_coord) {
		waiting_time = fmod(pickup_time_second_transfer + (p.nodes[node_destination_terminal].x_coord + p.nodes[node_destination_terminal].y_coord), pt_interval);
		if (waiting_time < -0.0001) { waiting_time = waiting_time + pt_interval; }
	}
	else {
		waiting_time = fmod(pickup_time_second_transfer + (pt_interval - p.nodes[node_destination_terminal].x_coord + p.nodes[node_destination_terminal].y_coord), pt_interval);
		if (waiting_time < -0.0001) { waiting_time = waiting_time + pt_interval; }
	}

	return travel_time + service_duration + waiting_time;
}


void update_load(struct problem &p, struct solution &s, int ref_first_node, int ref_last_node) {

	int ref = ref_first_node;
	int previous_ref = s.predecessor[ref];

	int vehicle = s.map_to_vehicle[ref];
	int start_depot = 2 * p.n_requests + vehicle;

	// No update needed for the start depot (if included):
	if (ref == start_depot) { previous_ref = ref; ref = s.successor[ref]; }

	while (ref != s.successor[ref_last_node]) {

		int request = s.map_to_request[ref];

		// Nodes where users are picked up:
		if (ref < p.n_requests || ref >= 3 * p.n_requests + 2 * p.n_vehicles) {
			s.load_mobile[ref] = s.load_mobile[previous_ref] + p.requests[request].users_mobile;
			s.load_wheelchair[ref] = s.load_wheelchair[previous_ref] + p.requests[request].users_wheelchair;
		}

		// Nodes where users are delivered:
		else if (ref < 2 * p.n_requests || ref >= 2 * p.n_requests + 2 * p.n_vehicles) {
			s.load_mobile[ref] = s.load_mobile[previous_ref] - p.requests[request].users_mobile;
			s.load_wheelchair[ref] = s.load_wheelchair[previous_ref] - p.requests[request].users_wheelchair;
		}
		
		previous_ref = ref;
		ref = s.successor[ref];
	}
}


void update_earliest_time(struct problem &p, struct solution &s, int begin_node) {

	int ref = begin_node;
	int node = s.map_to_node[ref];
	int previous_ref = s.predecessor[ref];
	int previous_node = s.map_to_node[previous_ref];

	// Update earliest start at begin depot (if included in the range):
	if (ref >= 2 * p.n_requests && ref < 2 * p.n_requests + p.n_vehicles) {
		s.earliest_time[ref] = p.nodes[node].lower_tw;
		previous_ref = ref;
		previous_node = node;
		ref = s.successor[ref];
		node = s.map_to_node[ref];
	}

	// Update earliest start at succeeding nodes including the end depot:
	while (true) {
		s.earliest_time[ref] = max(p.nodes[node].lower_tw, s.earliest_time[previous_ref] + p.nodes[previous_node].service_dur + p.matrix[previous_node * p.n_nodes + node]);
		if (ref >= 2 * p.n_requests && ref < 2 * p.n_requests + 2 * p.n_vehicles) { break; }
		previous_ref = ref;
		previous_node = node;
		ref = s.successor[ref];
		node = s.map_to_node[ref];
	}
}


void update_latest_time(struct problem &p, struct solution &s, int end_node) {

	int ref = end_node;
	int node = s.map_to_node[ref];
	int next_ref = s.successor[ref];
	int next_node = s.map_to_node[next_ref];

	// Update latest start at end depot (if included in the range):
	if (ref >= 2 * p.n_requests + p.n_vehicles && ref < 2 * p.n_requests + 2 * p.n_vehicles) {
		s.latest_time[ref] = p.nodes[node].upper_tw;
		next_ref = ref;
		next_node = node;
		ref = s.predecessor[ref];
		node = s.map_to_node[ref];
	}

	// Update latest start at preceding nodes including the begin depot:
	while (true) {
		s.latest_time[ref] = min(p.nodes[node].upper_tw, s.latest_time[next_ref] - p.nodes[node].service_dur - p.matrix[node * p.n_nodes + next_node]);
		if (ref >= 2 * p.n_requests && ref < 2 * p.n_requests + p.n_vehicles) { break; }
		next_ref = ref;
		next_node = node;
		ref = s.predecessor[ref];
		node = s.map_to_node[ref];
	}
}


bool check_feasibility(struct problem &p, struct solution &s) {

	// Feasibility checks related to requests:
	for (int request = 0; request < p.n_requests; request++) {

		int pickup = request;
		int delivery = request + p.n_requests;
		int ref_transfer1 = request + 2 * p.n_requests + 2 * p.n_vehicles;
		int ref_transfer2 = request + 3 * p.n_requests + 2 * p.n_vehicles;

		// Check that the request is served:
		if (s.map_to_vehicle[pickup] == -1 || s.map_to_vehicle[delivery] == -1) { printf("ERROR - Request is not served"); return false; }

		// Check time windows at pickup and delivery:
		if (s.service_start[pickup] < p.nodes[pickup].lower_tw - 0.001 || s.service_start[pickup] > p.nodes[pickup].upper_tw + 0.001) { printf("ERROR - Infeasible service time"); return false; }
		if (s.service_start[delivery] < p.nodes[delivery].lower_tw - 0.001 || s.service_start[delivery] > p.nodes[delivery].upper_tw + 0.001) { printf("ERROR - Infeasible service time"); return false; }

		// Check maximum user ride time:
		if (s.service_start[delivery] - (s.service_start[pickup] + p.nodes[pickup].service_dur) > p.requests[request].max_ride_time + 0.001) { printf("ERROR - Infeasible user ride time"); return false; }

		// If no transfers are made, check that pickup and delivery are served in the same route:
		if ((s.map_to_vehicle[ref_transfer1] == -1 || s.map_to_vehicle[ref_transfer2] == -1) && (s.map_to_vehicle[pickup] != s.map_to_vehicle[delivery])) { printf("ERROR - Infeasible routes assigned"); return false; }

		// If transfers are made, check the synchronization at the transfer locations:
		if (s.map_to_vehicle[pickup] != s.map_to_vehicle[delivery]) {
			if (s.service_start[ref_transfer2] < s.service_start[ref_transfer1] + compute_required_time_pt(p, s.map_to_node[ref_transfer1], s.map_to_node[ref_transfer2], s.service_start[ref_transfer2]) - 0.001) {
				printf("ERROR - Infeasible synchronization"); return false;
			}
		}
	}

	// Feasibility checks related to routes:
	for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {

		int begin_depot = 2 * p.n_requests + vehicle;
		int end_depot = begin_depot + p.n_vehicles;

		// Check maximum route duration:
		if (s.service_start[end_depot] - s.service_start[begin_depot] > p.vehicles[vehicle].max_route_dur + 0.001) { printf("ERROR - Infeasible route duration"); return false; }

		// Check travel times between successive nodes in a route:
		int ref = begin_depot;
		while (ref != end_depot) {
			if (s.service_start[ref] + p.nodes[s.map_to_node[ref]].service_dur + p.matrix[s.map_to_node[ref] * p.n_nodes + s.map_to_node[s.successor[ref]]] > s.service_start[s.successor[ref]] + 0.001) { printf("ERROR - Infeasible travel time"); return false; }
			ref = s.successor[ref];
		}

		// Check loads (with upgrading conditions):
		ref = s.successor[begin_depot];
		while (ref != end_depot) {

			// Check loads at the start depot:
			if (s.load_mobile[begin_depot] != 0 || s.load_wheelchair[begin_depot] != 0) {
				printf("ERROR - Infeasible load"); return false;
			}

			// Check loads at the end depot:
			if (s.load_mobile[end_depot] != 0 || s.load_wheelchair[end_depot] != 0) {
				printf("ERROR - Infeasible load"); return false;
			}

			// Nodes where users are picked up:
			if (ref < p.n_requests || ref >= 3 * p.n_requests + 2 * p.n_vehicles) {
				if ((s.load_mobile[s.predecessor[ref]] + p.requests[s.map_to_request[ref]].users_mobile != s.load_mobile[ref]) || (s.load_mobile[ref] > p.vehicles[s.map_to_vehicle[ref]].cap_mobile) ||
					(s.load_wheelchair[s.predecessor[ref]] + p.requests[s.map_to_request[ref]].users_wheelchair != s.load_wheelchair[ref]) || (s.load_wheelchair[ref] > p.vehicles[s.map_to_vehicle[ref]].cap_wheelchair)) {
					printf("ERROR - Infeasible load"); return false;
				}
			}

			// Nodes where users are delivered:		
			else if (ref < 2 * p.n_requests || ref >= 2 * p.n_requests + 2 * p.n_vehicles) {
				if ((s.load_mobile[s.predecessor[ref]] - p.requests[s.map_to_request[ref]].users_mobile != s.load_mobile[ref]) || (s.load_mobile[ref] > p.vehicles[s.map_to_vehicle[ref]].cap_mobile) ||
					(s.load_wheelchair[s.predecessor[ref]] - p.requests[s.map_to_request[ref]].users_wheelchair != s.load_wheelchair[ref]) || (s.load_wheelchair[ref] > p.vehicles[s.map_to_vehicle[ref]].cap_wheelchair)) {
					printf("ERROR - Infeasible load"); return false;
				}
			}

			ref = s.successor[ref];
		}
	}

	return true;
}


void write_output_file(struct problem &p, struct solution &s, int current_iteration, double computation_time, bool extensive_output) {

	// Compute output related to ride times:
	double actual_ride_time = 0, actual_ride_time_noPT = 0, actual_ride_time_PT = 0;
	double excess_ride_time = 0, excess_ride_time_noPT = 0, excess_ride_time_PT = 0;
	double detour_factor = 0, detour_factor_noPT = 0, detour_factor_PT = 0;
	double waiting_time_PT = 0;

	for (int request = 0; request < p.n_requests; request++) {

		int pickup = request;
		int delivery = pickup + p.n_requests;
		int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request;
		int ref_terminal2 = ref_terminal1 + p.n_requests;
		int terminal1 = s.map_to_node[ref_terminal1] - 2 * p.n_requests - p.n_vehicles;
		int terminal2 = s.map_to_node[ref_terminal2] - 2 * p.n_requests - p.n_vehicles;

		double request_direct_rt = p.matrix_nontightened[pickup * p.n_nodes + delivery];
		double request_actual_rt = s.service_start[delivery] - (s.service_start[pickup] + p.nodes[pickup].service_dur);
		double request_excess_rt = request_actual_rt - request_direct_rt;
		double request_detour_factor = request_actual_rt / request_direct_rt;

		actual_ride_time += request_actual_rt;
		excess_ride_time += request_excess_rt;
		detour_factor += request_detour_factor;
		
		if (s.map_to_vehicle[ref_terminal1] != -1 && transfers == 1) {
			actual_ride_time_PT += request_actual_rt;
			excess_ride_time_PT += request_excess_rt;
			detour_factor_PT += request_detour_factor;
			waiting_time_PT += (s.service_start[ref_terminal2] - (s.service_start[ref_terminal1] + p.nodes[s.map_to_node[ref_terminal1]].service_dur) - p.matrix_pt[terminal1 * p.n_terminals + terminal2]);
		}
		else {
			actual_ride_time_noPT += request_actual_rt;
			excess_ride_time_noPT += request_excess_rt;
			detour_factor_noPT += request_detour_factor;
		}
	}

	// Register users served by public transport:
	double count_pt = 0;
	vector<int> users_pt;
	if (transfers == 1) {
		for (int request = 0; request < p.n_requests; request++) {
			if (s.map_to_vehicle[2 * p.n_requests + 2 * p.n_vehicles + request] != -1) {
				count_pt = count_pt + 1;
				users_pt.push_back(request);
			}
		}
	}

	// Register usage of the public transport terminals:
	int *usage_terminals = new int[p.n_terminals];
	fill(usage_terminals, usage_terminals + p.n_terminals, 0);
	if (transfers == 1) {
		for (int request = 0; request < p.n_requests; request++) {
			if (s.map_to_vehicle[2 * p.n_requests + 2 * p.n_vehicles + request] != -1) {

				int ref_terminal1 = 2 * p.n_requests + 2 * p.n_vehicles + request;
				int ref_terminal2 = ref_terminal1 + p.n_requests;
				int terminal1 = s.map_to_node[ref_terminal1] - 2 * p.n_requests - p.n_vehicles;
				int terminal2 = s.map_to_node[ref_terminal2] - 2 * p.n_requests - p.n_vehicles;
				
				usage_terminals[terminal1] += 1;
				usage_terminals[terminal2] += 1;
			}
		}
	}


	ofstream output_file;

	// Write summary file:
	output_file.open(("results_summary.txt"), std::ios_base::app);
	output_file << data_darp << '\t'
				<< transfers << '\t'
				<< relative_mrt << '\t'
				<< mrt_factor << '\t'
				<< speed_factor_pt << '\t'
				<< pt_interval << '\t'
				<< iterations << '\t'
				<< iter_without_impr << '\t'
				<< max_removal_pct << '\t'
				<< max_deterior_pct << '\t'
				<< clustered_removal << '\t'
				<< idarp_adapted << '\t'
				<< without_operator << '\t'
				<< terminals_per_node << '\t'
				<< priority_non_pt << '\t'
				<< priority_pt << '\t'
				<< priority_long_dist << '\t'
				<< priority_short_dist << '\t'
				<< current_iteration << '\t'
				<< s.total_distance << '\t'
				<< computation_time << '\t'
				<< actual_ride_time / p.n_requests << '\t'
				<< actual_ride_time_PT / count_pt << '\t'
				<< actual_ride_time_noPT / (p.n_requests - count_pt) << '\t'
				<< excess_ride_time / p.n_requests << '\t'
				<< excess_ride_time_PT / count_pt << '\t'
				<< excess_ride_time_noPT / (p.n_requests - count_pt) << '\t'
				<< detour_factor / p.n_requests << '\t'
				<< detour_factor_PT / count_pt << '\t'
				<< detour_factor_noPT / (p.n_requests - count_pt) << '\t'
				<< waiting_time_PT / count_pt << '\t'
				<< count_pt / p.n_requests << '\t' 
				<< (double)accept_random_removal / max(1, current_iteration) << '\t'
				<< (double)accept_worst_removal / max(1, current_iteration) << '\t'
				<< (double)accept_related_removal / max(1, current_iteration) << '\t'
				<< (double)accept_route_removal / max(1, current_iteration) << '\t'
				<< (double)accept_random_order_insertion / max(1, current_iteration) << '\t'
				<< (double)accept_greedy_insertion / max(1, current_iteration) << '\t'
				<< (double)accept_two_regret_insertion / max(1, current_iteration) << '\t'
				<< (double)improve_random_removal / max(1, current_iteration) << '\t'
				<< (double)improve_worst_removal / max(1, current_iteration) << '\t'
				<< (double)improve_related_removal / max(1, current_iteration) << '\t'
				<< (double)improve_route_removal / max(1, current_iteration) << '\t'
				<< (double)improve_random_order_insertion / max(1, current_iteration) << '\t'
				<< (double)improve_greedy_insertion / max(1, current_iteration) << '\t'
				<< (double)improve_two_regret_insertion / max(1, current_iteration) << '\t';
	for (int m = 0; m < p.n_terminals; m++) { output_file << usage_terminals[m] << '\t'; }
	output_file << endl;
	output_file.close();

	delete[] usage_terminals;

	if (extensive_output == true) {

		// Write extensive output file:
		output_file.open("results_extensive.txt", std::ios_base::app);
		output_file << "Data DARP: " << data_darp << endl
					<< "Transfers: " << transfers << endl
					<< "Relative MRT: " << relative_mrt << endl
					<< "MRT factor: " << mrt_factor << endl
					<< "Speed factor PT : " << speed_factor_pt << endl
					<< "PT interval: " << pt_interval << endl
					<< "Iterations: " << iterations << endl
					<< "Iter without impr: " << iter_without_impr << endl
					<< "Max removal pct: " << max_removal_pct << endl
					<< "Max deterior pct: " << max_deterior_pct << endl
					<< "Clustered removal: " << clustered_removal << endl
					<< "IDARP adapted: " << idarp_adapted << endl
					<< "Without operator:  " << without_operator << endl
					<< "Terminals per node: " << terminals_per_node << endl
					<< "Priority non PT: " << priority_non_pt << endl
					<< "Priority PT: " << priority_pt << endl
					<< "Priority long: " << priority_long_dist << endl
					<< "Priority short: " << priority_short_dist << endl
					<< "Current iteration: " << current_iteration << endl
					<< "Total distance: " << s.total_distance << endl
					<< "Computation time: " << computation_time << endl
					<< "Average ride time: " << actual_ride_time / p.n_requests << endl
					<< "Pct PT users: " << count_pt / p.n_requests << endl;

		// List of users served by public transport:
		output_file << "Users PT: ";
		for (int user = 0; user < users_pt.size(); user++) {
			output_file << users_pt[user] << " ";
		}
		output_file << endl;

		// List of vehicle routes with load and time schedule:
		for (int vehicle = 0; vehicle < p.n_vehicles; vehicle++) {
			int ref = 2 * p.n_requests + vehicle;
			int end_depot = 2 * p.n_requests + p.n_vehicles + vehicle;
			output_file << endl << "- route " << vehicle << ": ";
			while (true) {
				output_file << ref << " (" << s.load_mobile[ref] << " " << s.load_wheelchair[ref] << " " << s.service_start[ref] << ") ";
				if (ref == end_depot) { break; }
				else { ref = s.successor[ref]; }
			}
			output_file << endl;
		}
		output_file << endl << endl;

		output_file.close();
	}
}