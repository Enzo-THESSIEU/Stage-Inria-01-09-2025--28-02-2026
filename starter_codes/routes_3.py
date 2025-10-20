def extract_route_3(vars_, nodes, N, K):
    v = vars_["v"]
    T_node = vars_["T_node"]
    B_node = vars_["B_node"]

    # Dictionary to store used arcs and their values
    used_arcs = {}

    # Extract all used arcs
    for (i, j) in v.keys():
        if v[(i, j)].X > 0.5:  # Only consider arcs that are used
            used_arcs[(i, j)] = [[i,j], i, nodes[nodes[:, 0] == i, 1][0], nodes[nodes[:, 0] == i, 2][0] if len(nodes[nodes[:, 0] == i, 2]) > 0 else "", T_node[j].X, B_node[j].X]
    print(used_arcs)

    used_arcs_vehicle_1 = []
    used_arcs_vehicle_2 = []
    first_trips = []

    keys = list(used_arcs.keys())
    for key in keys:
        if key[0] == 0:  # Vehicle 1 depot
            first_trips.append(used_arcs.pop(key))

    iter = 0

    print("First trips:", first_trips)
    print(len(first_trips))
    # assert len(first_trips) == 2, "Error: Number of first trips does not match number of vehicles"
    if len(first_trips) == 2:
        used_arcs_vehicle_1.append(first_trips[0])
        used_arcs_vehicle_2.append(first_trips[1])
        # while len(used_arcs) > 0:
        while iter < 1000:
            iter += 1
            for key in keys:
                if key[0] == used_arcs_vehicle_1[-1][0][1]:  # Continue vehicle 1
                    used_arcs_vehicle_1.append(used_arcs.pop(key))
                if key[0] == used_arcs_vehicle_2[-1][0][1]:  # Continue vehicle 2
                    used_arcs_vehicle_2.append(used_arcs.pop(key))
            # print(used_arcs_vehicle_1)
            # print(used_arcs_vehicle_2)

    elif len(first_trips) == 1:
        used_arcs_vehicle_1.append(first_trips[0])
        # while len(used_arcs) > 0:
        while iter < 1000:
            iter += 1
            for key in keys:
                if key[0] == used_arcs_vehicle_1[-1][0][1]:  # Continue vehicle 1
                    used_arcs_vehicle_1.append(used_arcs.pop(key))
            # print(len(used_arcs))
            # print(iter)



    return used_arcs_vehicle_1, used_arcs_vehicle_2

def extract_route_4(vars_, nodes, N, K, zeroDepot_node, endDepot_node):
    v = vars_["v"]
    T_node = vars_["T_node"]
    B_node = vars_["B_node"]

    used_arcs = {}

    for (i, j) in v.keys():
        if v[(i, j)].X > 0.5:
            # Handle tuple node IDs
            base_i = i[0] if isinstance(i, tuple) else i
            base_j = j[0] if isinstance(j, tuple) else j
            print("base_i :", base_i)

            # Find node type and request info from nodes array
            node_type = nodes[nodes[:, 0] == base_i, 1][0]
            print("node_type :",node_type)
            req_info  = nodes[nodes[:, 0] == base_i, 2][0] if len(nodes[nodes[:, 0] == base_i, 2]) > 0 else ""
            print("req_info :", req_info)


            used_arcs[(i, j)] = [
                [i, j],
                base_i,
                node_type,
                req_info,
                T_node[j].X,
                B_node[j].X
            ]

    routes = {k: [] for k in K}

    for k in K:
        current_node = zeroDepot_node
        visited = set()
        while True:
            found_next = False
            for (i, j), info in used_arcs.items():
                if i == current_node and (i, j) not in visited:
                    routes[k].append(info)
                    visited.add((i, j))
                    current_node = j
                    found_next = True
                    if j == endDepot_node:
                        break
            if not found_next or current_node == endDepot_node:
                break

    return routes
