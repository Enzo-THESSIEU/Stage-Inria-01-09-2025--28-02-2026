// define the header for small instance generator code
#ifndef SMALL_INSTANCE_GENERATOR_H
#define SMALL_INSTANCE_GENERATOR_H

// Include necessary libraries
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>

// Include the different classes used

// --------------------
// Request_Generator
// --------------------
class Request_Generator{
    private:
        int first_node_location;
        int second_node_location;
        double pi = 3.141592653589793;

        std::vector<double> R;      // [R1, R2, R3]
        std::vector<double> r;      // [r_1, r_2]
        std::vector<double> theta;  // [theta_1, theta_2]
        double max_time = 1380;     // in minutes (23 hours)
        double min_time = 480;      // in minutes (8 hours)
        double latest_arrival_time; // Latest arrival time for the request

            // ---- Random r ----
        double sample_r(std::mt19937& gen, double R) {
            static std::uniform_real_distribution<double> U(0.0, 1.0);
            double u = U(gen);
            return -R * std::log(1.0 - u * (1.0 - std::exp(-1.0)));
        }

        // ---- Random theta ----
        double sample_theta(std::mt19937& gen){        
            std::uniform_real_distribution<double> theta_dist(0.0, 2.0 * pi);
            return theta_dist(gen);
        }

        double sample_arrival_time(std::mt19937& gen, double min_time, double max_time) {
            std::uniform_real_distribution<double> time_dist(min_time, max_time);
            return time_dist(gen);
        }

    public:
        Request_Generator();
        int get_first_node_location() const;
        int get_second_node_location() const;
        double get_r1() const;
        double get_r2() const;
        double get_theta1() const;
        double get_theta2() const;
        double get_latest_arrival_time() const;
        
        void print() const;



};


// --------------------
// Timetabled_Train_Generator
// --------------------
class Timetabled_Train_Generator {
    public:
        struct Train {
            int id;
            std::string stations;
            double departure;
            double arrival;
        };

        std::vector<std::vector<Train>> timetable;

        Timetabled_Train_Generator();

        std::vector<Train> get_all_trains() const;
        void print_timetable() const;

    private:

        // PT arcs: (1<->2), (1<->3), (2<->3)
        std::vector<double> pt_travel_time        = {60.0, 75.0, 45.0};
        std::vector<double> pt_earliest_departure = {390.0, 390.0, 390.0};
        std::vector<double> pt_latest_departure   = {1290.0, 1290.0, 1290.0};
        std::vector<double> pt_frequency          = {90.0, 75.0, 60.0};
        std::vector<std::string> arc_name         = {"(1<->2)", "(1<->3)", "(2<->3)"};

};


// --------------------
// All_Requests_Generator
// --------------------
class All_Requests_Generator {
    public:
        explicit All_Requests_Generator(int n_requests);
        void print_requests() const;

    private:
        struct request {
            int request_id;
            int pickup_node_id;
            int dropoff_node_id;
            int pickup_node_hub;
            int dropoff_node_hub;
            double earliest_pickup_time;
            double latest_arrival_time;
            double service_time;
            double x_coord_pickup;
            double y_coord_pickup;
            double x_coord_dropoff;
            double y_coord_dropoff;
        };

        int n_requests;
        std::vector<Request_Generator> request_generated;
        std::vector<request> all_requests;
        Timetabled_Train_Generator train_generator;

};

#endif