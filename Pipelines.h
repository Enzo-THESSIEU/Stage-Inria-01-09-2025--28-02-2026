#pragma once

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

// Structures
#include "structures.h"

// Pipelines
void import_results_pipeline();
void import_most_correlated_requests_pipiline();
void update_TW_pipeline();
void build_TW_score_for_subset();
void build_request_db_pipeline();
std::string print_n_tw_score_requests(int instance_id, int subset_number, int size_of_subset, int size_of_set, int offset);

instance_info get_instance_info_pipeline(int instance_id);





