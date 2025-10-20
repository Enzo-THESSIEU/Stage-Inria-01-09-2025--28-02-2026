##### Imports #####
import time as t
import gurobipy as gb

### Original Model and Arc Elimination ###
from Problem_definition_1 import build_idarp_model_Original, build_idarp_model_Enhanced_1
from routes_1 import extract_routes
import Graph_builder_1 as g1

### All Other Enhancements ###
from Problem_definition_2 import build_idarp_model_Enhanced_2, build_idarp_model_Enhanced_3, build_idarp_model_Enhanced_4
from routes_2 import extract_route_2
import Graph_Builder_2 as g2 
from Cluster_Heuristic import subtour_callback

from Problem_Definition_3 import build_idarp_model_Enhanced_5_EV
import Graph_Builder_3 as g3
from routes_3 import extract_route_3, extract_route_4

from Problem_Definition_4 import build_idarp_model_Enhanced_6_Timetabled_Departures
import Graph_Builder_4 as g4

Time_Limit = 10#* 60 * 60

tij = g3.data['tij_original']
P = g3.data['P']
D = g3.data['D']
n = len(P)
pair_pi_di = g3.data['pair_pi_di']

def excess_user_ride_time(route):
    service_times = {}
    total_excess = 0
    for node in route:
        service_time_node = node[4]
        if node in P + D:
            service_times[node] = service_time_node
    for key in pair_pi_di:
        total_ride_time = service_times[pair_pi_di[key]] - service_times[pair_pi_di[key]]
        excess_user_ride_time_val = total_ride_time - tij[key, pair_pi_di[key]] 
        total_excess += excess_user_ride_time_val
    return total_excess

def Build_route_Original(Time_Limit):
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Original(Time_Limit)
    m.optimize()
    end_cpu = t.process_time()
    end = t.perf_counter()
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_routes(
        vars_,
        nodes=g1.data['nodes'],
        N=g1.data['N'],
        K=g1.data['K'],
        zeroDepot=g1.data['zeroDepot'],
        endDepot=g1.data['endDepot']
    )
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_1(Time_Limit):
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_1(Time_Limit)
    m.optimize()
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_routes(
        vars_,
        nodes=g1.data['nodes'],
        N=g1.data['N'],
        K=g1.data['K'],
        zeroDepot=g1.data['zeroDepot'],
        endDepot=g1.data['endDepot']
    )
    end_cpu = t.process_time()
    end = t.perf_counter()
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_2(Time_Limit): 
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_2(Time_Limit)

    m.optimize()
    end_cpu = t.process_time()
    end = t.perf_counter()
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_route_2(vars_, nodes=g2.data['nodes'], N=g2.data['N'], K=g2.data['K'])
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_3(Time_Limit):
    start_cpu = t.process_time() 
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_3(Time_Limit)

    m.optimize(subtour_callback)
    end_cpu = t.process_time()
    end = t.perf_counter()
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_route_2(vars_, nodes=g2.data['nodes'], N=g2.data['N'], K=g2.data['K'])
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_4(Time_Limit): 
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_4(Time_Limit)
    m.optimize(subtour_callback)
    end_cpu = t.process_time()
    end = t.perf_counter()

    # Extract and print routes
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_route_2(vars_, nodes=g2.data['nodes'], N=g2.data['N'], K=g2.data['K'])
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_5_EV(Time_Limit): 
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_5_EV(Time_Limit)
    m.optimize(subtour_callback)
    end_cpu = t.process_time()
    end = t.perf_counter()

    # Extract and print routes
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_route_3(vars_, nodes=g3.data['nodes'], N=g3.data['N'], K=g3.data['K'])
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles

def Build_route_Enhanced_6_Timetabled_Departures(Time_Limit):
    start_cpu = t.process_time()
    start = t.perf_counter()
    m, vars_ = build_idarp_model_Enhanced_6_Timetabled_Departures(Time_Limit)
    m.optimize(subtour_callback)
    end_cpu = t.process_time()
    end = t.perf_counter()

    # Extract and print routes
    used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_route_4(vars_, nodes=g4.data['nodes'], N=g4.data['N'], K=g4.data['K'], zeroDepot_node = g4.data['zeroDepot_node'], endDepot_node = g4.data['endDepot_node'])
    excess_time_vehicles = excess_user_ride_time(used_arcs_vehicle_1) - excess_user_ride_time(used_arcs_vehicle_2)
    time_elapsed_cpu = start_cpu - end_cpu
    elapsed = end - start
    return m.ObjVal, m.ObjBound, used_arcs_vehicle_1, used_arcs_vehicle_2, -time_elapsed_cpu, elapsed, excess_time_vehicles


# === Run all models and store results ===
results = {
    "Original":       Build_route_Original(Time_Limit),
    "Enhanced_1":     Build_route_Enhanced_1(Time_Limit),
    "Enhanced_2":     Build_route_Enhanced_2(Time_Limit),
    "Enhanced_3":     Build_route_Enhanced_3(Time_Limit),
    "Enhanced_4":     Build_route_Enhanced_4(Time_Limit),
    "Enhanced_5_EV":  Build_route_Enhanced_5_EV(Time_Limit),
    "Enhanced_6_Timetabled_Departures": Build_route_Enhanced_6_Timetabled_Departures(Time_Limit)
}

print(results)

# === Pretty print results ===
# for name, (objval, bound, v1_route, v2_route, cpu_time, elapsed) in results.items():
#     print("---------------------------------------------------------------")
#     print(f"{name} Model Objective Value:", objval)
#     print(f"{name} Model Lower Bound:", bound)
#     print(f"{name} Vehicle 1 Route:", v1_route)
#     print(f"{name} Vehicle 2 Route:", v2_route)
#     print(f"{name} CPU Solve Time (s):", cpu_time)
#     print(f"{name} Elapsed Solve Time (s):", elapsed)







