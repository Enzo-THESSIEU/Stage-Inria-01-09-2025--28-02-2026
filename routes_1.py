def extract_routes(vars_, nodes, N, K, zeroDepot, endDepot):
    x       = vars_["x"]        # keyed by (k,i,j)
    B_node  = vars_["B_node"]
    B_veh   = vars_["B_veh"]

    def arc_record(i, j, time_at_j):
        node_type = nodes[nodes[:, 0] == i, 1][0]
        req_cell  = nodes[nodes[:, 0] == i, 2]
        request   = req_cell[0] if len(req_cell) > 0 else ""
        return [[i, j], i, node_type, request, time_at_j]

    # Build the route for a given vehicle k by following x[k,*,*]
    def build_route_for_vehicle(k):
        route = []

        # find first hop from depot
        start_js = [j for j in N if (k, zeroDepot, j) in x and x[k, zeroDepot, j].X > 0.5]
        if not start_js:
            return route  # empty route if vehicle unused
        j0 = start_js[0]

        # depot departure record (match routes_2 shape; time = B_veh[k])
        route.append(arc_record(zeroDepot, j0, B_veh[k].X))

        # follow until endDepot
        cur = j0
        visited_guard = set()
        while cur != endDepot:
            visited_guard.add(cur)
            next_js = [j for j in N if (k, cur, j) in x and x[k, cur, j].X > 0.5]
            if not next_js:
                break
            j1 = next_js[0]
            # time at j1 uses B_node[j1]
            route.append(arc_record(cur, j1, B_node[j1].X))
            cur = j1
            # safety: break cycles if any
            if cur in visited_guard and cur != endDepot:
                break

        return route

    # Assume exactly two vehicles for parity with routes_2 usage
    v1 = build_route_for_vehicle(K[0]) if len(K) > 0 else []
    v2 = build_route_for_vehicle(K[1]) if len(K) > 1 else []

    return v1, v2
