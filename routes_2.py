class DARPRouteExtractor:
    def __init__(self, m, vars_, sets, params, **options):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params

        # Unpack all boolean flags, with defaults
        self.duplicate_transfers = options.get("duplicate_transfers", True)
        self.arc_elimination = options.get("arc_elimination", True)
        self.variable_substitution = options.get("variable_substitution", True)
        self.subtour_elimination = options.get("subtour_elimination", True)
        self.transfer_node_strengthening = options.get("transfer_node_strengthening", True)
        self.ev_constraints = options.get("ev_constraints", False)
        self.timetabled_departures = options.get("timetabled_departures", False)
        self.use_imjn = options.get("use_imjn", False)
        self.MoPS = options.get("MoPS", False)


        

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i


    def safe_node_info(self, i):
        """Safely fetch node info (type, request id) for node i."""
        i_base = self.base(i)
        nodes = self.sets['nodes']
        try:
            row = nodes[nodes[:, 0] == i_base][0]
            node_type = row[1]
            req_id = row[2] if len(row) > 2 else ""
        except Exception:
            node_type, req_id = "unknown", ""
        return node_type, req_id


    def extract_vehicle_route_final(self):
        """
        Universal route extractor — works with:
        - variable_substitution: True (uses v) or False (uses x)
        - use_imjn: True (artificial node tuples) or False (flat node IDs)
        Orders request routes by service time (T_node).
        """
        nodes, N, P, D, C, F, R, K, zeroDepot_node, endDepot_node, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        y, z, T_node = self.vars_["y"], self.vars_['z'], self.vars_["T_node"]

        # === VEHICLE ARCS ===
        used_arcs = {}

        if self.variable_substitution:
            v = self.vars_["v"]
            for (i, j) in v.keys():
                if v[(i, j)].X > 0.5:
                    node_type, req_id = self.safe_node_info(i)
                    used_arcs[(i, j)] = [
                        [i, j],
                        i,
                        node_type,
                        req_id,
                        f"Service time at node {i}",
                        T_node[i].X if (i in T_node) else None,
                    ]
        else:
            x = self.vars_["x"]
            for (k, i, j) in x.keys():
                if x[(k, i, j)].X > 0.5:
                    node_type, req_id = self.safe_node_info(i)
                    used_arcs[(i, j)] = [
                        [i, j],
                        i,
                        node_type,
                        req_id,
                        f"Service time at node {i}",
                        T_node[i].X if (i in T_node) else None,
                    ]

        # --- Detect first trips (from depot) ---
        used_arcs_vehicle_1, used_arcs_vehicle_2 = [], []
        first_trips = [used_arcs.pop(k) for k in list(used_arcs.keys()) if self.base(k[0]) == 0]

        if first_trips:
            used_arcs_vehicle_1.append(first_trips[0])
        if len(first_trips) > 1:
            used_arcs_vehicle_2.append(first_trips[1])

        # --- Stitch vehicle paths ---
        stuck_counter = 0
        while used_arcs and stuck_counter < 100:
            progress = False
            for key in list(used_arcs.keys()):
                i, j = key
                if used_arcs_vehicle_1 and i == used_arcs_vehicle_1[-1][0][1]:
                    used_arcs_vehicle_1.append(used_arcs.pop(key))
                    progress = True
                elif used_arcs_vehicle_2 and i == used_arcs_vehicle_2[-1][0][1]:
                    used_arcs_vehicle_2.append(used_arcs.pop(key))
                    progress = True
            if not progress:
                stuck_counter += 1
                if stuck_counter == 3:
                    print("⚠️ No progress in vehicle stitching — breaking loop.")
                    break
        return used_arcs_vehicle_1, used_arcs_vehicle_2, used_arcs

    def extract_request_route_final(self):

        nodes, N, P, D, C, F, R, K, zeroDepot_node, endDepot_node, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        y, z, T_node = self.vars_["y"], self.vars_['z'], self.vars_["T_node"]

        # === REQUEST ARCS ===
        used_arcs_y = {
            (int(r), i, j): [int(r), [i, j]]
            for (r, i, j) in y.keys()
            if y[(r, i, j)].X > 0.5
        }

        # Group arcs per request
        req_routes = {r: [] for r in set(r for (r, _, _) in used_arcs_y.keys())}

        # Identify earliest node for each request
        for (r, i, j), val in list(used_arcs_y.items()):
            if not req_routes[r] or (
                T_node[i].X if i in T_node else float("inf")
            ) < (T_node[req_routes[r][0][1][0]].X if req_routes[r] else float("inf")):
                req_routes[r].insert(0, val)
                used_arcs_y.pop((r, i, j), None)

        # --- Stitch sequentially ---
        iter_limit = 200
        for _ in range(iter_limit):
            if not used_arcs_y:
                break
            progress = False
            for key in list(used_arcs_y.keys()):
                r, i, j = key
                if req_routes[r]:
                    last_j = req_routes[r][-1][1][1]
                    if i == last_j:
                        req_routes[r].append(used_arcs_y.pop(key))
                        progress = True
            if not progress:
                break

        # Append leftover arcs (if disjoint)
        if used_arcs_y:
            print(f"⚠️ Unconnected arcs left for some requests: {len(used_arcs_y)}")
            for key in used_arcs_y.keys():
                r, i, j = key
                req_routes.setdefault(r, []).append(used_arcs_y[key])

        # === SORT REQUEST ROUTES BY SERVICE TIME ===
        for r in req_routes:
            req_routes[r].sort(
                key=lambda a: T_node[a[1][1]].X if a[1][1] in T_node else float("inf")
            )

        return (
            req_routes.get(1, []),
            req_routes.get(2, []),
            req_routes.get(3, []),
            req_routes.get(4, []),
        )
    
    def test_passenger_transfer_balance(self, route1, route2):
        nodes, N, P, D, C, F, R, K = (self.sets[k] for k in ["nodes","N","P","D","C","F","R","K"])
        zeroDepot_node, endDepot_node, A = self.sets["zeroDepot"], self.sets["endDepot"], [tuple(a) for a in self.sets["A"]]
        y, z, T_node = self.vars_["y"], self.vars_["z"], self.vars_["T_node"]
        fi_r, Departures = self.params["fi_r"], self.params["Departures"]

        sum_balance = {}

        nodes_visited = []

        for route in (route1, route2):
            for step in route:
                if isinstance(step, list) and len(step) > 1:
                    node = step[1]  # e.g. (11,1)
                    nodes_visited.append(node)

        for r in R:
            for node in nodes_visited:
                sum_yij = sum(int(y[r, node, j].X) for j in N if (node, j) in A)
                sum_yji = sum(int(y[r, j, node].X) for j in N if (j, node) in A)
                sum_zij = 0
                sum_zji = 0

                if node in C:
                    for j in C:
                        if (node, j) in A and (self.base(node), self.base(j)) in Departures:
                            for d in Departures[(self.base(node), self.base(j))]:
                                sum_zij += int(z[d, node, j].X)
                        if (j, node) in A and (self.base(j), self.base(node)) in Departures:
                            for d in Departures[(self.base(j), self.base(node))]:
                                sum_zji += int(z[d, j, node].X)

                lhs = sum_yij + sum_zij - sum_yji - sum_zji
                rhs = fi_r[(r, node)]
                sum_balance[(r, node)] = {
                    "y_out": sum_yij,
                    "z_out": sum_zij,
                    "y_in": sum_yji,
                    "z_in": sum_zji,
                    "lhs": lhs,
                    "rhs": rhs,
                    "diff": lhs - rhs
                }

        return sum_balance

    def extract_PT_route_final(self):
        """
        Extracts all used Public Transport (PT) arcs from the model.
        Handles both timetabled and non-timetabled cases robustly.
        """
        nodes = self.sets["nodes"]
        C = self.sets["C"]
        z = self.vars_["z"]
        T_node = self.vars_["T_node"]
        used_PT_arcs = []

        print(f"PT extraction mode: {'timetabled' if self.timetabled_departures else 'simple'}")
        print(f"→ z variable keys example: {list(z.keys())[:5]}")

        if self.timetabled_departures:
            Departures = self.params.get("Departures", {})
            print(f"→ Departures keys: {list(Departures.keys())[:5]}")

            for i in C:
                for j in C:
                    key = (self.base(i), self.base(j))
                    if key not in Departures:
                        continue

                    for d in Departures[key]:
                        if (d, i, j) in z:
                            if z[d, i, j].X > 1e-6:
                                used_PT_arcs.append([
                                    (i, j),
                                    f"Departure {d}",
                                    f"T({i})={T_node[i].X:.2f}",
                                    f"T({j})={T_node[j].X:.2f}",
                                    f"z={z[d,i,j].X:.3f}"
                                ])
                        elif (i, j) in z:  # fallback if not 3D
                            if z[i, j].X > 1e-6:
                                used_PT_arcs.append([
                                    (i, j),
                                    f"T({i})={T_node[i].X:.2f}",
                                    f"T({j})={T_node[j].X:.2f}",
                                    f"z={z[i,j].X:.3f}"
                                ])
        else:
            for i in C:
                for j in C:
                    if (i, j) in z and z[i, j].X > 1e-6:
                        used_PT_arcs.append([
                            (i, j),
                            f"T({i})={T_node[i].X:.2f}",
                            f"T({j})={T_node[j].X:.2f}",
                            f"z={z[i,j].X:.3f}"
                        ])

        print(f"🚌 Extracted {len(used_PT_arcs)} PT arcs.")
        return used_PT_arcs




            # === Return ordered routes ===
            
