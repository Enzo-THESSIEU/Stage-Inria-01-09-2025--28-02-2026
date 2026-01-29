#pragma once

// Necessary Libraries
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


// Structures
#include "structures.h"


class table_creators{
    public:
    void create_results_table_if_missing(const sql_params& sql_param);
    void create_request_table_if_missing(const sql_params& sql_param);
    void delete_table(const sql_params& sql_param, const std::string& table);
};


class table_modifiers{
    public:
    void add_column_to_table(const sql_params& sql_param,
                         const std::string& table,
                         const std::string& newColumn,
                         const std::string& columnType);
};


class info_getters{
    public:
    std::vector<int> get_instances(sql_params sql_param);
    instance_info get_instance_info(int instance_id, sql_params sql_param);
    std::pair<std::vector<node_info>, std::vector<node_info>> link_pickups_dropoffs(int instance_id, instance_info instance, sql_params sql_param);
    node_info get_node_info(int instance_id, int node_id, sql_params sql_param);
    std::vector<node_info> get_nodes_by_type(int instance_id, const std::string& node_type, int limit, sql_params sql_param);


    double avg_distance(int iteration, sql_params sql_param);
    double avg_ride_time(int iteration, sql_params sql_param);
    double avg_detour_factor(int iteration, sql_params sql_param);
    
    std::vector<db_request_info> get_n_requests(int n, int offset, int instance_id, sql_params sql_param);
};


class input_to_db{
    public: 
    void import_results_summary_tsv(const std::string& path, int instance_id, sql_params& sql_param);
    void input_constrained_TW_into_db(std::vector<node_info> pickup_nodes, std::vector<node_info> dropoff_nodes, int instance_id, sql_params& sql_param);
    int insert_trains_into_db(const std::vector<Train>& all_trains, sql_params& sql_param);
    void insert_correlated_set(const sql_params& sql_param, std::vector<int> correlated_set, int instance_id, int request_id);
    void input_to_db_tw_scores(const sql_params& sql_param, double tw_score, int instance_id, int request_id);
    void input_to_db_request_info(const sql_params& sql_param, db_request_info req);
};