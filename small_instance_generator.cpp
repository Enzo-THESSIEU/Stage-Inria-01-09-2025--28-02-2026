#include "small_instance_generator.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>


Request_Generator::Request_Generator() {

        // ---- Fixed radii ----
        R = {20.0, 17.0, 19.0};

        // ---- Random engine ----
        static std::random_device rd;
        static std::mt19937 gen(rd());

        // ---- Choose node locations ----
        std::vector<int> locations = {1, 2, 3};
        std::shuffle(locations.begin(), locations.end(), gen);

        first_node_location  = locations[0];
        second_node_location = locations[1];

        r = {
            sample_r(gen, R[first_node_location - 1]),
            sample_r(gen, R[second_node_location - 1])
        };

        theta = {
            sample_theta(gen),
            sample_theta(gen)
        };

        latest_arrival_time = sample_arrival_time(gen, min_time, max_time);
    }

int Request_Generator::get_first_node_location() const {
    return first_node_location;
}

int Request_Generator::get_second_node_location() const {
    return second_node_location;
}

double Request_Generator::get_r1() const {
    return r[0];
}

double Request_Generator::get_r2() const {
    return r[1];
}

double Request_Generator::get_theta1() const {
    return theta[0];
}

double Request_Generator::get_theta2() const {
    return theta[1];
}

double Request_Generator::get_latest_arrival_time() const {
    return latest_arrival_time;
}

void Request_Generator::print() const {
    std::cout << "First node: " << first_node_location << "\n";
    std::cout << "Second node: " << second_node_location << "\n";

    std::cout << "r: ";
    for (double x : r) std::cout << x << " ";
    std::cout << "\n";

    std::cout << "theta: ";
    for (double x : theta) std::cout << x << " ";
    std::cout << "\n";

    std::cout << "Latest arrival time: " << latest_arrival_time << "\n";
}



Timetabled_Train_Generator::Timetabled_Train_Generator() {

    int n_arcs = pt_travel_time.size();
    timetable.resize(n_arcs);
    int train_id = 0;

    for (int i = 0; i < n_arcs; i++){


        for (double departure = pt_earliest_departure[i]; departure <= pt_latest_departure[i]; departure += pt_frequency[i]){
            double arrival = departure + pt_travel_time[i];
            std::string stations = arc_name[i];
            timetable[i].push_back({train_id, stations, departure, arrival});
            train_id++;
        }
    }      
}

void Timetabled_Train_Generator::print_timetable() const {

    for (size_t arc = 0; arc < timetable.size(); ++arc) {
        std::cout << "Arc " << arc_name[arc] << "\n";

        for (const auto& tr : timetable[arc]) {
            std::cout << "  Train Id: " << tr.id << "\n";
            std::cout << "  Stations: " << tr.stations << "\n";

            std::cout << "  Departure (min): " << tr.departure << "\n";
            std::cout << "  Departure: " << (int)(tr.departure / 60)
                    << "h" << (int)tr.departure % 60 << "\n";

            std::cout << "  Arrival (min): " << tr.arrival << "\n";
            std::cout << "  Arrival: " << (int)(tr.arrival / 60)
                    << "h" << (int)tr.arrival % 60 << "\n";

            std::cout << "---------------------\n";
        }
    }
}

std::vector<Timetabled_Train_Generator::Train> Timetabled_Train_Generator::get_all_trains() const {
    std::vector<Train> out;
    for (const auto& arcVec : timetable) {
        out.insert(out.end(), arcVec.begin(), arcVec.end());
    }
    return out;
}

All_Requests_Generator::All_Requests_Generator(int n_requests){
    for (int request_id = 0; request_id < n_requests; request_id++) {

            // Generate request parameters
            request_generated.emplace_back();
            const Request_Generator& gen = request_generated.back();

            int first_node_location  = gen.get_first_node_location();
            int second_node_location = gen.get_second_node_location();

            double r_1 = gen.get_r1();
            double r_2 = gen.get_r2();
            double theta_1 = gen.get_theta1();
            double theta_2 = gen.get_theta2();
            double latest_arrival_time = gen.get_latest_arrival_time();

            request formatted_request;

            formatted_request.request_id = request_id;
            formatted_request.pickup_node_id = request_id;
            formatted_request.dropoff_node_id = request_id + n_requests;
            formatted_request.pickup_node_hub = first_node_location;
            formatted_request.dropoff_node_hub = second_node_location;

            formatted_request.latest_arrival_time = latest_arrival_time;
            formatted_request.earliest_pickup_time = 0;
                // latest_arrival_time
                // - 1.5 * (
                //     t[first_node_location][second_node_location]
                // + t[request_id][first_node_location]
                // + t[second_node_location][request_id + n_requests]
                // );

            formatted_request.service_time = 1.0;

            formatted_request.x_coord_pickup  = r_1 * std::cos(theta_1);
            formatted_request.y_coord_pickup  = r_1 * std::sin(theta_1);
            formatted_request.x_coord_dropoff = r_2 * std::cos(theta_2);
            formatted_request.y_coord_dropoff = r_2 * std::sin(theta_2);

            all_requests.push_back(formatted_request);
    }
}

void All_Requests_Generator::print_requests() const {
    std::cout << "Generated " << all_requests.size() << " requests\n\n";

    for (const auto& r : all_requests) {
        std::cout << "Request " << r.request_id << "\n";
        std::cout << "  Pickup node ID:   " << r.pickup_node_id << "\n";
        std::cout << "  Dropoff node ID:  " << r.dropoff_node_id << "\n";
        std::cout << "  Pickup hub:       " << r.pickup_node_hub << "\n";
        std::cout << "  Dropoff hub:      " << r.dropoff_node_hub << "\n";
        std::cout << "  Earliest pickup:  " << r.earliest_pickup_time << "\n";
        std::cout << "  Latest arrival:   " << r.latest_arrival_time << "\n";
        std::cout << "  Service time:     " << r.service_time << "\n";
        std::cout << "  Pickup coords:    ("
                << r.x_coord_pickup << ", "
                << r.y_coord_pickup << ")\n";
        std::cout << "  Dropoff coords:   ("
                << r.x_coord_dropoff << ", "
                << r.y_coord_dropoff << ")\n";
        std::cout << "-----------------------------\n";
    }
}
    
