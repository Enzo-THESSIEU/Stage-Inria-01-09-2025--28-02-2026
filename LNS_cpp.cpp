#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <variant>
#include <optional>
#include <algorithm>
#include <limits>
#include <tuple>

// Type aliases for readability
using Node = std::variant<int, std::tuple<int, int>>;  // can represent an int or (i,j)
using Matrix = std::map<std::pair<int, int>, double>;  // for tij, etc.
using Dict = std::map<int, double>;

class InitialSolution {
public:
    // === Attributes ===
    std::map<std::string, std::vector<int>> sets;
    std::map<std::string, std::variant<Matrix, Dict, std::map<int, int>, int>> params;

    bool duplicate_transfers;
    bool arc_elimination;
    bool variable_substitution;
    bool subtour_elimination;
    bool transfer_node_strengthening;
    bool ev_constraints;
    bool timetabled_departures;
    bool use_imjn;
    bool MoPS;

    int n_veh;
    std::vector<std::vector<int>> vehicle_route;
    std::vector<std::vector<int>> vehicle_current_requests;
    std::vector<int> loc;
    std::vector<double> load;

    // === Constructor ===
    InitialSolution(
        const std::map<std::string, std::vector<int>>& sets_,
        const std::map<std::string, std::variant<Matrix, Dict, std::map<int, int>, int>>& params_,
        const std::map<std::string, bool>& options = {}
    ) : sets(sets_), params(params_)
    {
        duplicate_transfers = get_option(options, "duplicate_transfers", true);
        arc_elimination = get_option(options, "arc_elimination", true);
        variable_substitution = get_option(options, "variable_substitution", true);
        subtour_elimination = get_option(options, "subtour_elimination", true);
        transfer_node_strengthening = get_option(options, "transfer_node_strengthening", true);
        ev_constraints = get_option(options, "ev_constraints", false);
        timetabled_departures = get_option(options, "timetabled_departures", false);
        use_imjn = get_option(options, "use_imjn", false);
        MoPS = get_option(options, "MoPS", false);

        n_veh = static_cast<int>(sets.at("K").size());

        int zeroDepot = std::get<int>(params.at("zeroDepot_node"));

        vehicle_route = std::vector<std::vector<int>>(n_veh, std::vector<int>{zeroDepot});
        vehicle_current_requests = std::vector<std::vector<int>>(n_veh);
        loc = std::vector<int>(n_veh, zeroDepot);
        load = std::vector<double>(n_veh, 0.0);
    }

    // === Helper to extract base of node ===
    int base(const Node& i) const {
        if (std::holds_alternative<std::tuple<int, int>>(i))
            return std::get<0>(std::get<std::tuple<int, int>>(i));
        else
            return std::get<int>(i);
    }

    // === Set base extraction ===
    std::set<int> base_set(const std::set<Node>& node_list) const {
        std::set<int> res;
        for (const auto& i : node_list)
            res.insert(base(i));
        return res;
    }

    // === Sort by earliest time window ===
    std::vector<int> order_by_earliest_time_window(const std::vector<int>& unsorted_set) {
        const auto& ei = std::get<Dict>(params.at("ei"));
        std::vector<int> sorted_set = unsorted_set;
        std::sort(sorted_set.begin(), sorted_set.end(),
            [&](int a, int b) { return ei.at(base(a)) < ei.at(base(b)); });
        return sorted_set;
    }

    // === Find closest dropoff ===
    std::optional<int> closest_dropoff(int v) {
        double min_time = std::numeric_limits<double>::infinity();
        std::optional<int> closest_drop;

        const auto& tij = std::get<Matrix>(params.at("tij"));
        const auto& pair_pi_di = std::get<std::map<int, int>>(params.at("pair_pi_di"));

        for (int req : vehicle_current_requests[v]) {
            int drop_node = pair_pi_di.at(req);
            double current_time = tij.at({ base(loc[v]), base(drop_node) });
            if (current_time < min_time) {
                min_time = current_time;
                closest_drop = drop_node;
            }
        }
        return closest_drop;
    }

    // === Find closest vehicle to a pickup ===
    int closest_vehicle(int p) {
        const auto& tij = std::get<Matrix>(params.at("tij"));
        double dist_closest = std::numeric_limits<double>::infinity();
        int closest = -1;
        for (int v = 0; v < n_veh; ++v) {
            double dist = tij.at({ base(loc[v]), base(p) });
            if (dist < dist_closest) {
                dist_closest = dist;
                closest = v;
            }
        }
        return closest;
    }

    // === Algorithm 2: Closest Request Greedy Heuristic ===
    void closest_request_greedy_heuristic() {
        auto& tij = std::get<Matrix>(params.at("tij"));
        auto unserved_pickups = order_by_earliest_time_window(sets.at("P"));

        while (!unserved_pickups.empty()) {
            int p = unserved_pickups.back();
            unserved_pickups.pop_back();

            int v = closest_vehicle(p);
            auto closest_drop = closest_dropoff(v);

            if (!closest_drop || tij[{ base(loc[v]), base(p) }] <
                                 tij[{ base(loc[v]), base(*closest_drop) }]) {
                loc[v] = p;
                vehicle_route[v].push_back(p);
            } else {
                loc[v] = *closest_drop;
                vehicle_route[v].push_back(*closest_drop);
                unserved_pickups.push_back(p);
            }
        }
    }

    // === Algorithm 3: Earliest Time Window Greedy Heuristic ===
    void earliest_request_greedy_heuristic() {
        auto& tij = std::get<Matrix>(params.at("tij"));
        auto unserved_P = order_by_earliest_time_window(sets.at("P"));
        auto unserved_D = order_by_earliest_time_window(sets.at("D"));

        std::vector<int> unserved_P_and_D = unserved_P;
        unserved_P_and_D.insert(unserved_P_and_D.end(), unserved_D.begin(), unserved_D.end());

        const auto& pair_pi_di = std::get<std::map<int, int>>(params.at("pair_pi_di"));

        while (!unserved_P_and_D.empty()) {
            int next_node = unserved_P_and_D.back();
            unserved_P_and_D.pop_back();

            if (std::find(unserved_D.begin(), unserved_D.end(), next_node) != unserved_D.end()) {
                // It's a dropoff node
                int pickup_node = -1;
                for (const auto& kv : pair_pi_di)
                    if (kv.second == next_node)
                        pickup_node = kv.first;

                int vehicle_index = -1;
                for (int v = 0; v < n_veh; ++v)
                    if (std::find(vehicle_route[v].begin(), vehicle_route[v].end(), pickup_node) != vehicle_route[v].end())
                        vehicle_index = v;

                if (vehicle_index >= 0)
                    vehicle_route[vehicle_index].push_back(next_node);
            } else {
                // It's a pickup
                int next_vehicle = closest_vehicle(next_node);
                vehicle_route[next_vehicle].push_back(next_node);
                loc[next_vehicle] = next_node;
            }
        }
    }

private:
    // === Helper to extract options with default ===
    bool get_option(const std::map<std::string, bool>& options, const std::string& key, bool default_val) {
        auto it = options.find(key);
        return (it != options.end()) ? it->second : default_val;
    }
};
