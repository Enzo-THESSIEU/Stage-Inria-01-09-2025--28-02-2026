class DARPDebuggingFunctions:
    """Utility class for debugging and analyzing DARP model outputs."""

    def __init__(self, m, vars_, sets, params, **options_1):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params
        self.base = lambda i: i[0] if isinstance(i, tuple) else i

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

    def extract_variable_values(self):
        """Extract and return all main variable values."""
        z_values_all = {k: int(v.X) for k, v in self.vars_['z'].items()}
        z_values_1 = {k: int(v.X) for k, v in self.vars_['z'].items() if abs(v.X) > 1e-6}

        y_values_r4 = {k: int(v.X) for k, v in self.vars_['y'].items() if k[0] == 4}

        a_values = {}
        if self.use_imjn:
            a_values = {k: int(v.X) for k, v in self.vars_['a'].items()}

        v_values_1 = {}
        sum_v_values = {}
        if self.variable_substitution:
            v_values_1 = {
                k: self.params['tij'][self.base(k[0]),self.base(k[1])]
                for k, v in self.vars_['v'].items()
                if abs(v.X) > 1e-6
            }
            sum_v_values = sum(v_values_1.values())


        x_values_1 = {k: int(v.X) for k, v in self.vars_['x'].items() if abs(v.X) > 1e-6}

        return z_values_all, z_values_1, y_values_r4, a_values, v_values_1, sum_v_values, x_values_1

    def check_transfer_balance(self, extractor):
        """Check PT transfer balance using extractor."""
        sum_balance = extractor.test_passenger_transfer_balance(*extractor.extract_vehicle_route_final())
        wrong_sum_balance = {
            k: v for k, v in sum_balance.items() if abs(v['diff']) > 1e-6
        }
        return wrong_sum_balance
    
    def check_flow_conservation(self):
        P, D, C, N, K, A = self.sets['P'], self.sets['D'], self.sets['C'], self.sets['N'], self.sets['K'], self.sets['A']
        x= self.vars_['x']
        flow_conservation = {}
        problems = []
        for k in K:
            for i in P+D+C:
                inflow = sum(x[k, j, i].X for j in N if (j,i) in x.keys())
                outflow = sum(x[k, i, j].X for j in N if (i,j) in x.keys())
                flow_conservation[f"{k}", f"{i}"] = {
                    "inflow": inflow,
                    "outflow": outflow,
                    "flow conservation": inflow - outflow
                }
                if inflow - outflow != 0:
                    problems.append([k, i])
        return flow_conservation, problems


    def compute_constraint_balance(self):
        """Compute passenger flow balance for each request and node."""
        R, C, N, A = self.sets['R'], self.sets['C'], self.sets['N'], self.sets['A']
        P, D = self.sets['P'], self.sets["D"]
        Departures = self.params['Departures']
        fi_r = self.params['fi_r']
        z = self.vars_['z']

        Constraint_val = {}
        problems = []
        for r in R:
            for i in P+D+C:
                # Passenger movements
                sum_pickup = sum(
                    self.vars_['y'][r, i, j].X for j in N if (i, j) in A and (r, i, j) in self.vars_['y']
                )
                sum_dropoff = sum(
                    self.vars_['y'][r, j, i].X for j in N if (j, i) in A and (r, j, i) in self.vars_['y']
                )

                # PT transfers
                z_sum_pickup, z_sum_dropoff = 0, 0
                if self.timetabled_departures:
                    for j in N:
                        if (self.base(i), self.base(j)) in Departures:
                            for d in Departures[(self.base(i), self.base(j))]:
                                if (d, i, j) in z.keys():
                                    z_sum_pickup += z[(d, i, j)].X
                        if (self.base(j), self.base(i)) in Departures:
                            for d in Departures[(self.base(j), self.base(i))]:
                                if (d, j, i) in z:
                                    z_sum_dropoff += z[(d, j, i)].X

                else:
                    for j in C:
                        if (i,j) in A:
                            if (i,j) in z.keys():
                                z_sum_pickup += z[(i, j)].X
                        if (j,i) in A:
                            if (j,i) in z.keys():
                                z_sum_dropoff += z[(j, i)].X

                lhs = sum_pickup + z_sum_pickup - sum_dropoff - z_sum_dropoff
                rhs = fi_r.get((r, i), 0)
                if lhs - rhs != 0:
                    problems.append([(int(r), i), sum_pickup, sum_dropoff, z_sum_pickup, z_sum_dropoff, lhs, rhs, lhs - rhs])

                Constraint_val[(int(r), i)] = {
                    'y_out': sum_pickup,
                    'y_in': sum_dropoff,
                    'z_out': z_sum_pickup,
                    'z_in': z_sum_dropoff,
                    'lhs': lhs,
                    'rhs': rhs,
                    'diff': lhs - rhs,
                }
        return Constraint_val, problems
    
    def print_z_variable_summary(self, model_name):
        """
        Print a summary of z variables, check for missing PT arcs,
        and print final objective info.
        """
        vars_ = self.vars_
        C = self.sets["C"]
        Departures = self.params["Departures"]

        Total_z_variables = len(vars_['z'])
        first_ten_elements = list(vars_['z'].keys())[:10]

        count_existing = 0
        count_missing = 0

        if self.timetabled_departures:
            for i in C:
                for j in C:
                    if (self.base(i), self.base(j)) in Departures:
                        for d in Departures[(self.base(i), self.base(j))]:
                            if (d, i, j) in vars_['z']:
                                count_existing += 1
                            else:
                                count_missing += 1
        
        else: 
            for i in C:
                for j in C:
                    if (i,j) in self.sets['A']:
                        if (i, j) in vars_['z']:
                            count_existing += 1
                        else:
                            count_missing += 1

        return Total_z_variables, first_ten_elements, count_existing, count_missing

class DARPRouteDebugging:
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

    def extract_vehicle_route_final_ev(self):
        """
        Universal route extractor.
        Returns enriched route lists for Veh1 and Veh2:
        [
        [ arc=(i,j), node=i, type, req, T_i, SOC%, charge_time, pct_charged ]
        ...
        ]
        """

        # === Shortcuts ===
        zeroDepot = self.sets["zeroDepot"]
        endDepot = self.sets["endDepot"]
        K = self.sets["K"]
        F = self.sets["F"]

        T_node = self.vars_["T_node"]
        B_node = self.vars_["B_node"]
        E_node = self.vars_["E_node"]

        C_bat = self.params["C_bat"]
        alpha = self.params["alpha"]

        # -----------------------------------------
        # Helper
        # -----------------------------------------
        def soc(i):
            """Battery % at node i."""
            if i in B_node:
                # return round(B_node[i].X / C_bat * 100, 2)
                return round(B_node[i].X / C_bat * 100, 3)
            return None

        def charge_info(i, node_type):
            """Return (charging_time, percent_charged)."""
            if node_type != "charging station":
                return 0.0, 0.0
            if i in E_node:
                ct = float(E_node[i].X)
                pct = round(alpha * ct / C_bat * 100, 1)
                return ct, pct
            return 0.0, 0.0

        # -----------------------------------------
        # CASE 1 – variable substitution v[i,j]
        # -----------------------------------------
        if self.variable_substitution:
            v = self.vars_["v"]

            # collect used arcs
            arcs = [(i, j) for (i, j) in v.keys() if v[i, j].X > 0.5]

            # separate by vehicle using depot start
            veh_routes = {0: [], 1: []}
            unassigned = arcs.copy()

            # detect arcs starting from depot
            starts = [arc for arc in arcs if self.base(arc[0]) == 0]

            if len(starts) >= 1:
                veh_routes[0].append(starts[0])
                unassigned.remove(starts[0])
            if len(starts) >= 2:
                veh_routes[1].append(starts[1])
                unassigned.remove(starts[1])

            # follow paths
            for k in [0, 1]:
                route = veh_routes[k]
                stuck = 0
                while unassigned and stuck < 50 and route:
                    last = route[-1][1]
                    found = False
                    for arc in list(unassigned):
                        if arc[0] == last:
                            route.append(arc)
                            unassigned.remove(arc)
                            found = True
                            break
                    if not found:
                        stuck += 1
                veh_routes[k] = route

            # enrich
            enriched_routes = {0: [], 1: []}
            for k in [0, 1]:
                out = []
                for (i, j) in veh_routes[k]:
                    node_type, req = self.safe_node_info(i)
                    T = T_node[i].X if i in T_node else None
                    SoC = soc(i)
                    ct, pct = charge_info(i, node_type)
                    out.append([[i, j], i, node_type, req, T, SoC, ct, pct])
                enriched_routes[k] = out

            return enriched_routes[0], enriched_routes[1]

        # -----------------------------------------
        # CASE 2 – classical x[k,i,j]
        # -----------------------------------------
        x = self.vars_["x"]
        veh_routes = {k: [] for k in K}

        for k in K:
            used = [(i, j) for (kk, i, j) in x.keys() if kk == k and x[kk, i, j].X > 0.5]

            if not used:
                continue

            # find start arc
            try:
                first_arc = next((i, j) for (i, j) in used if i == zeroDepot)
            except StopIteration:
                continue

            route = [first_arc]
            used.remove(first_arc)

            # follow chain
            stuck = 0
            while stuck < 200:
                last_j = route[-1][1]
                if last_j == endDepot:
                    break
                found = False
                for arc in used:
                    if arc[0] == last_j:
                        route.append(arc)
                        used.remove(arc)
                        found = True
                        break
                if not found:
                    stuck += 1
                else:
                    stuck = 0

            veh_routes[k] = route

        # -----------------------------------------
        # enrich & return in Veh1, Veh2 order
        # -----------------------------------------
        enriched_routes = {}
        for k in K:
            out = []
            for (i, j) in veh_routes[k]:
                node_type, req = self.safe_node_info(i)
                T = T_node[i].X if i in T_node else None
                SoC = soc(i)
                ct, pct = charge_info(i, node_type)
                out.append([[i, j], i, node_type, req, T, SoC, ct, pct])
            enriched_routes[k] = out

        return enriched_routes.get(0, []), enriched_routes.get(1, [])

    def extract_vehicle_route_final(self):
        """
        Universal route extractor WITHOUT any EV constraints.
        Returns enriched route lists for Veh1 and Veh2:
        [
            [ arc=(i,j), node=i, type, req, T_i ]
        ]
        """

        # === Shortcuts ===
        zeroDepot = self.sets["zeroDepot"]
        endDepot = self.sets["endDepot"]
        K = self.sets["K"]

        T_node = self.vars_.get("T_node", {})

        # -----------------------------------------
        # Helper: extract T_i safely
        # -----------------------------------------
        def get_time(i):
            if i in T_node:
                try:
                    return float(T_node[i].X)
                except:
                    return None
            return None

        # -----------------------------------------
        # CASE 1 – variable substitution v[i,j]
        # -----------------------------------------
        if self.variable_substitution:
            v = self.vars_.get("v", {})

            # collect used arcs
            arcs = [(i, j) for (i, j) in v.keys() if v[i, j].X > 0.5]

            veh_routes = {0: [], 1: []}
            unassigned = arcs.copy()

            # detect arcs starting at depot
            starts = [arc for arc in arcs if self.base(arc[0]) == 0]

            if len(starts) >= 1:
                veh_routes[0].append(starts[0])
                unassigned.remove(starts[0])

            if len(starts) >= 2:
                veh_routes[1].append(starts[1])
                unassigned.remove(starts[1])

            # follow chains
            for k in [0, 1]:
                route = veh_routes[k]
                stuck = 0
                while unassigned and stuck < 50 and route:
                    last = route[-1][1]
                    found = False
                    for arc in list(unassigned):
                        if arc[0] == last:
                            route.append(arc)
                            unassigned.remove(arc)
                            found = True
                            break
                    if not found:
                        stuck += 1
                veh_routes[k] = route

            # enrich
            enriched = {0: [], 1: []}
            for k in [0, 1]:
                out = []
                for (i, j) in veh_routes[k]:
                    node_type, req = self.safe_node_info(i)
                    T = get_time(i)
                    out.append([
                        [i, j],
                        i,
                        node_type,
                        req,
                        T
                    ])
                enriched[k] = out

            return enriched[0], enriched[1]

        # -----------------------------------------
        # CASE 2 – classical x[k,i,j]
        # -----------------------------------------
        x = self.vars_.get("x", {})
        veh_routes = {k: [] for k in K}

        for k in K:
            used = [(i, j) for (kk, i, j) in x.keys() if kk == k and x[kk, i, j].X > 0.5]

            if not used:
                continue

            try:
                first_arc = next((i, j) for (i, j) in used if i == zeroDepot)
            except StopIteration:
                continue

            route = [first_arc]
            used.remove(first_arc)

            # follow arc chain
            stuck = 0
            while stuck < 200:
                last_j = route[-1][1]
                if last_j == endDepot:
                    break

                found = False
                for arc in used:
                    if arc[0] == last_j:
                        route.append(arc)
                        used.remove(arc)
                        found = True
                        break

                if not found:
                    stuck += 1
                else:
                    stuck = 0

            veh_routes[k] = route

        # -----------------------------------------
        # Enrich nodes
        # -----------------------------------------
        enriched = {}
        for k in K:
            out = []
            for (i, j) in veh_routes[k]:
                node_type, req = self.safe_node_info(i)
                T = get_time(i)
                out.append([
                    [i, j],
                    i,
                    node_type,
                    req,
                    T
                ])
            enriched[k] = out

        return enriched.get(0, []), enriched.get(1, [])


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
                sum_yij = sum(int(y[r, node, j].X) for j in N if (r, node, j) in y.keys())
                sum_yji = sum(int(y[r, j, node].X) for j in N if (r, j, node) in y.kesy())
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

    def extract_PT_route_final(self):
        """
        Extracts all used Public Transport (PT) arcs from the model.
        Handles both timetabled and non-timetabled cases robustly.
        """
        nodes = self.sets["nodes"]
        R = self.sets['R']
        C = self.sets["C"]
        z = self.vars_["z"]
        T_node = self.vars_["T_node"]
        Cr = self.params['Cr']
        used_PT_arcs = []

        # print(f"PT extraction mode: {'timetabled' if self.timetabled_departures else 'simple'}")
        # print(f"→ z variable keys example: {list(z.keys())[:5]}")

        if self.timetabled_departures:
            Departures = self.params.get("Departures", {})
            # print(f"→ Departures keys: {list(Departures.keys())[:5]}")
            for r in R:
                for i in C:
                    for j in C:
                        key = (self.base(i), self.base(j))
                        if key not in Departures:
                            continue

                        for d in Departures[key]:
                            if (d, r, i, j) in z:
                                if z[d, r, i, j].X > 1e-6:
                                    used_PT_arcs.append([
                                        (i, j),
                                        f"Departure {d}",
                                        f"T({i})={T_node[i].X:.2f}",
                                        f"T({j})={T_node[j].X:.2f}",
                                        f"z={z[d,r,i,j].X:.3f}"
                                    ])

        else:
            if Cr:
                for r in R:
                    for i in Cr[r]:
                        for j in Cr[r]:
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
            
    def ev_constraints_issue(self):
        # === Shortcuts ===
        zeroDepot = self.sets["zeroDepot"]
        endDepot = self.sets["endDepot"]
        K = self.sets["K"]
        F = self.sets["F"]

        T_node = self.vars_["T_node"]
        B_node = self.vars_["B_node"]
        E_node = self.vars_["E_node"]

        C_bat = self.params["C_bat"]
        alpha = self.params["alpha"]
        beta  = self.params['beta']
        tij   = self.params['tij']

        ev_constraint_list = []

        if self.variable_substitution:
            # ------------------------------------------
            # CASE 1: Using v[i,j]
            # ------------------------------------------
            v = self.vars_["v"]

            for (i, j), var in v.items():
                if var.X > 0.5:
                    ev_constraint = []
                    ev_constraint.append(f"ARC USED: {(i,j)}")

                    # Battery values
                    Bi = round(B_node[i].X, 3)
                    Bj = round(B_node[j].X, 3)
                    ev_constraint.append(f"  Battery at {i}: {Bi}")
                    ev_constraint.append(f"  Battery at {j}: {Bj}")

                    # Travel consumption
                    t_ij = tij[(self.base(i), self.base(j))]
                    theoretical = round(beta * t_ij, 3)
                    actual = round(Bj - Bi, 3)
                    constraint_val = round((Bj - Bi + beta * t_ij), 3)

                    ev_constraint.append(f"  Theoretical consumption: {theoretical}")
                    ev_constraint.append(f"  Actual ΔB: {actual}")
                    ev_constraint.append(f"  Constraint: {constraint_val}")

                    ev_constraint_list.append(ev_constraint)

        else:
            # ------------------------------------------
            # CASE 2: Using x[k,i,j]
            # ------------------------------------------
            x = self.vars_["x"]

            # Loop over all arcs (k,i,j)
            for (k, i, j), var in x.items():
                if var.X > 0.5:
                    ev_constraint = []
                    ev_constraint.append(f"ARC USED: vehicle {k}, {(i,j)}")

                    # Battery values
                    Bi = round(B_node[i].X, 3)
                    Bj = round(B_node[j].X, 3)
                    ev_constraint.append(f"  Battery at {i}: {Bi}")
                    ev_constraint.append(f"  Battery at {j}: {Bj}")

                    # Travel consumption
                    t_ij = tij[(self.base(i), self.base(j))]
                    theoretical = round(beta * t_ij, 3)
                    actual = round(Bj - Bi, 3)
                    constraint_val = round((Bj - Bi + beta * t_ij), 3)

                    ev_constraint.append(f"  Theoretical consumption: {theoretical}")
                    ev_constraint.append(f"  Actual ΔB: {actual}")
                    ev_constraint.append(f"  Constraint: {constraint_val}")

                    ev_constraint_list.append(ev_constraint)

        return ev_constraint_list

