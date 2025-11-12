class DARPRouteExtractor:
    def __init__(self, m, vars_, sets, params, **options_1):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params
        
        options = options_1['options']

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


    def extract_vehicle_routes(self):
        """Extract vehicle routes from x[k,i,j]."""
        x = self.vars_.get("x", {})
        K = self.sets["K"]
        A = self.sets["A"]
        zeroDepot = self.sets["zeroDepot"]
        endDepot = self.sets["endDepot"]

        routes = {}
        arcs_use_k = [[(k2, i, j) for (k2, i, j) in x.keys() if x[k2, i, j].X > 0.5]]
        for k in K:
            arcs_used = [(i, j) for (k2, i, j) in x.keys() if k2 == k and x[k2, i, j].X > 0.5]
            if not arcs_used:
                routes[k] = []
                continue

            start_arc = next(((i, j) for (i, j) in arcs_used if i == zeroDepot), None)
            if not start_arc:
                print(f"[Warning] Vehicle {k} has no outgoing arc from depot.")
                routes[k] = []
                continue

            route = [start_arc]
            current = start_arc[1]
            arcs_used.remove(start_arc)

            while current != endDepot and arcs_used:
                next_arc = next(((i, j) for (i, j) in arcs_used if i == current), None)
                if not next_arc:
                    print(f"⚠️ Vehicle {k} route interrupted at {current}")
                    break
                route.append(next_arc)
                current = next_arc[1]
                arcs_used.remove(next_arc)

            routes[k] = route
        return routes, arcs_use_k


    def extract_request_routes(self):
        """Extract request routes from y[r,i,j]."""
        y = self.vars_.get("y", {})
        T_node = self.vars_.get("T_node", {})
        R = self.sets["R"]

        routes = {r: [] for r in R}
        for (r, i, j) in y.keys():
            if y[r, i, j].X == 1:
                routes[r].append((i, j))

        for r in R:
            routes[r].sort(key=lambda arc: T_node[arc[0]].X if arc[0] in T_node else 0)
        return routes

    
    def test_passenger_transfer_balance(self, route1, route2):
        nodes, N, P, D, C, F, R, K = (self.sets[k] for k in ["nodes","N","P","D","C","F","R","K"])
        zeroDepot_node, endDepot_node, A = self.sets["zeroDepot"], self.sets["endDepot"], [tuple(a) for a in self.sets["A"]]
        y, z, T_node = self.vars_["y"], self.vars_["z"], self.vars_["T_node"]
        fi_r, Departures = self.params["fi_r"], self.params["Departures"]
        Cr = self.params['Cr']

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

                if self.timetabled_departures:

                    if node in C:
                        for j in C:
                            if (node, j) in A and (self.base(node), self.base(j)) in Departures:
                                for d in Departures[(self.base(node), self.base(j))]:
                                    sum_zij += int(z[d, node, j].X)
                            if (j, node) in A and (self.base(j), self.base(node)) in Departures:
                                for d in Departures[(self.base(j), self.base(node))]:
                                    sum_zji += int(z[d, j, node].X)

                else:
                    if node in Cr[r]:
                        for j in Cr[r]:
                            if (node, j) in A:
                                sum_zij += int(z[node, j].X)
                            if (j, node) in A:
                                sum_zji += int(z[j, node].X)

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

    def extract_PT_routes(self):
        """Extract all used PT arcs z[d,i,j]."""
        z = self.vars_.get("z", {})
        TD, TA = self.params.get("TD", {}), self.params.get("TA", {})

        used_PT = []
        for key in z.keys():
            if len(key) == 3:
                d, i, j = key
                if z[d, i, j].X > 0.5:
                    dep = TD.get((self.base(i), self.base(j), d))
                    arr = TA.get((self.base(i), self.base(j), d))
                    used_PT.append(((i, j), d, dep, arr))
        return used_PT
    
    def summarize(self):
        veh_routes, arcs_use_k = self.extract_vehicle_routes()
        req_routes = self.extract_request_routes()
        pt_routes = self.extract_PT_routes()

        print("\n=== VEHICLE ROUTES ===")
        for k, route in veh_routes.items():
            if route:
                print(f"Vehicle {k}: " + " → ".join([f"{i}->{j}" for (i, j) in route]))
            else:
                print(f"Vehicle {k}: [no route]")
        print(arcs_use_k)

            

        print("\n=== REQUEST ROUTES ===")
        for r, arcs in req_routes.items():
            if arcs:
                print(f"Request {r}: " + " → ".join([f"{i}->{j}" for (i, j) in arcs]))
            else:
                print(f"Request {r}: [no assigned path]")

        print("\n=== PUBLIC TRANSPORT ARCS USED ===")
        if not pt_routes:
            print("No PT arcs used.")
        else:
            for ((i, j), d, dep, arr) in pt_routes:
                print(f"{i}->{j} (dep={dep}, arr={arr})")

        return veh_routes, req_routes, pt_routes






            # === Return ordered routes ===
            
