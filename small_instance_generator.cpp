#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

class Request_Generator {

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
        theta = { theta_dist(gen)};
        return theta_dist(gen);
    }

    double sample_arrival_time(std::mt19937& gen, double min_time, double max_time) {
        std::uniform_real_distribution<double> time_dist(min_time, max_time);
        return time_dist(gen);
    }

public:
    Request_Generator() {

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

        latest_arrival_time = {
            sample_arrival_time(gen, min_time, max_time)
        };
    }

    int get_first_node_location() const {
        return first_node_location;
    }

    int get_second_node_location() const {
        return second_node_location;
    }

    double get_r1() const {
        return r[0];
    }

    double get_r2() const {
        return r[1];
    }

    double get_theta1() const {
        return theta[0];
    }

    double get_theta2() const {
        return theta[1];
    }

    double get_latest_arrival_time() const {
        return latest_arrival_time;
    }

    void print() const {
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
};

class Timetabled_Train_Generator{
    private:

        struct Train {
            int id;
            double departure;
            double arrival;
        };

        // PT arcs: (1<->2), (1<->3), (2<->3)
        std::vector<double> pt_travel_time       = {60.0, 75.0, 45.0};
        std::vector<double> pt_earliest_departure = {390.0, 390.0, 390.0};
        std::vector<double> pt_latest_departure   = {1290.0, 1290.0, 1290.0};
        std::vector<double> pt_frequency          = {90.0, 75.0, 60.0};

        std::vector<std::vector<Train>> timetable;

    public:
    Timetabled_Train_Generator() {

        int n_arcs = pt_travel_time.size();
        timetable.resize(n_arcs);

        for (int i = 0; i < n_arcs; i++){
            int train_id = 0;

            for (double departure = pt_earliest_departure[i]; departure <= pt_latest_departure[i]; departure += pt_frequency[i]){
                double arrival = departure + pt_travel_time[i];
                timetable[i].push_back({train_id, departure, arrival});
                train_id++;
            }
        }      
    }
};

class All_Requests_Generator {
    private:
        struct request{
            int request_id; 
            int pickup_node_id; 
            int dropoff_node_id; 
            int pickup_node_hub; 
            int dropoff_node_hub; 
            double earliest_pickup_time; 
            double latest_arrival_time; 
            // int capacity_demand; 
            double service_time; 
            double x_coord_pickup;
            double y_coord_pickup;
            double x_coord_dropoff;
            double y_coord_dropoff;
        };
        // struct generated_request{
        //     int first_node_location;
        //     int second_node_location;
        //     double r_1;
        //     double r_2;
        //     double theta_1;
        //     double theta_2;
        //     double latest_arrival_time;
        // };

        int n_requests;
        std::vector<Request_Generator> request_generated;
        std::vector<request> all_requests;
        Timetabled_Train_Generator train_generator;

    public:
        All_Requests_Generator(int n_requests){
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



        void print_requests() const {
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
    
};

int main() {

    int n_requests = 6;
    int n_vehicles = 3;     // not yet used, but declared for clarity
    int n_pt_stations = 3;  // implicit in generators

    std::cout << "Generating small test instance\n";
    std::cout << "Requests: " << n_requests << "\n";
    std::cout << "Vehicles: " << n_vehicles << "\n";
    std::cout << "PT stations: " << n_pt_stations << "\n\n";

    All_Requests_Generator generator(n_requests);

    generator.print_requests();

    return 0;
};

