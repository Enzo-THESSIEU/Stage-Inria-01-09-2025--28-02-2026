##### Imports #####
import gurobipy as gb

class DARPConstraintBuilder:
    def __init__(self, m, vars_, sets, params, duplicate_transfers=True, arc_elimination=True, variable_substitution=True, subtour_elimination=True, transfer_node_strengthening=True,
                 ev_constraints=False, timetabled_departures=False, use_imjn=False, MoPS = False):
        self.duplicate_transfers = duplicate_transfers
        self.arc_elimination = arc_elimination
        self.variable_substitution = variable_substitution
        self.subtour_elimination = subtour_elimination
        self.transfer_node_strengthening = transfer_node_strengthening
        self.ev_constraints = ev_constraints
        self.timetabled_departures = timetabled_departures
        self.use_imjn = use_imjn
        self.MoPS = MoPS
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i

    def define_variables(self):

        """
        Define decision variables for the IDARP model, adapting to:
        - imjn artificial nodes (i,m)
        - timetabled departures (with z[d,i,j])
        - EV constraints (battery and charging)
        - variable substitution (link x and v)
        """

        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        Cr = self.params['Cr']
        Departures = self.params.get("Departures", None)

        ### --- Original variables ---
        # Vehicle routing arc variables
        if self.duplicate_transfers:
            transfer_arcs = {(0,0)}
            for r in R:
                for i in Cr[r]:
                    for j in Cr[r]:
                        transfer_arcs.add((i,j))
            transfer_arcs.remove((0,0))
        else:
            transfer_arcs = {(i,j) for i in C for j in C}

        DAR_arcs = {(i,j) for (i,j) in A if (i,j) not in transfer_arcs}
        request_arcs = {(i,j) for (i,j) in DAR_arcs if not (i == zeroDepot or j == endDepot)}

        x = self.m.addVars(K, DAR_arcs, vtype=gb.GRB.BINARY, name="x")   # vehicle arcs

        # Passenger flow
        y = self.m.addVars(R, request_arcs, vtype=gb.GRB.CONTINUOUS, lb=0, ub=1.0, name="y")
        # Node timing variables
        T_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")
        T_veh  = self.m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")
        D_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

        ### --- Developments ---
        # Optional variables
        v, a, B_node, E_node = None, None, None, None

        # Variable Substitution
        if self.variable_substitution:
            v = self.m.addVars(DAR_arcs, vtype=gb.GRB.BINARY, name="v")      # DAR arcs

        # Transfer / fixed-route arcs
        z = {}
        # count_z_keys_created = 0
        if self.timetabled_departures: # and Departures is not None:
            # Timetabled departures case
            for (i, j) in A:
                if i in C and j in C and (self.base(i), self.base(j)) in Departures:
                    for d in Departures[(self.base(i), self.base(j))]:
                        for r in R:
                            z[(d, r, i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{d},{i},{j}]")
                        # count_z_keys_created += 1
                        # print("z_keys_vounter : ", count_z_keys_created)
                        # print("Actual keys in z", len(z.keys()))
                        # print("\n")
                # Artificial node marker (used in some strengthened formulations)

        else:
            # Transfer Original Modeling
            for r in R:
                for i in Cr[r]:
                    for j in Cr[r]:
                        z[(i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")


        if self.use_imjn:
            a = self.m.addVars(N, vtype=gb.GRB.BINARY, name="a")


        # EV-specific variables
        if self.ev_constraints:
            B_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")  # battery SOC
            E_node = self.m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")  # recharge time

        self.vars_ = {"x": x, "v": v, "y": y, "z": z, "T_node": T_node, "T_veh": T_veh, "D_node": D_node, "a": a, "B_node": B_node, "E_node": E_node}
        return self.vars_

    # def define_variables(self):

    #         """
    #         Define decision variables for the IDARP model, adapting to:
    #         - imjn artificial nodes (i,m)
    #         - timetabled departures (with z[d,i,j])
    #         - EV constraints (battery and charging)
    #         - variable substitution (link x and v)
    #         """

    #         nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]

    #         Departures = self.params.get("Departures", None)

    #         ### --- Original variables ---
    #         # Vehicle routing arc variables
    #         x = self.m.addVars(K, A, vtype=gb.GRB.BINARY, name="x")   # vehicle arcs
    #         # Passenger flow
    #         y = self.m.addVars(R, A, vtype=gb.GRB.CONTINUOUS, lb=0, name="y")
    #         # Node timing variables
    #         T_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")
    #         T_veh  = self.m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")
    #         D_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

    #         ### --- Developments ---
    #         # Optional variables
    #         v, a, B_node, E_node = None, None, None, None
    #         z = {}

    #         # Variable Substitution
    #         if self.variable_substitution:
    #             v = self.m.addVars(A, vtype=gb.GRB.BINARY, name="v")      # DAR arcs

    #         # Transfer / fixed-route arcs
    #         z = {}
    #         if self.timetabled_departures and Departures is not None:
    #             # Timetabled departures case
    #             for (i, j) in A:
    #                 if i in C and j in C and (self.base(i), self.base(j)) in Departures:
    #                     for d in Departures[(self.base(i), self.base(j))]:
    #                         z[(d, i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{d},{i},{j}]")
    #                 # Artificial node marker (used in some strengthened formulations)

    #         if self.use_imjn:
    #             a = self.m.addVars(N, vtype=gb.GRB.BINARY, name="a")

    #         elif C:
    #             # Transfer Original Modeling
    #             for i in C:
    #                 for j in C:
    #                     z[(i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")

    #         # EV-specific variables
    #         if self.ev_constraints:
    #             B_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")  # battery SOC
    #             E_node = self.m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")  # recharge time

    #         return {
    #             "x": x, "v": v, "y": y, "z": z,
    #             "T_node": T_node, "T_veh": T_veh, "D_node": D_node,
    #             "a": a, "B_node": B_node, "E_node": E_node
    #         }


    # === Objective ===
    def set_objective(self, w = [1,0,0,0]):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        cij, tij, di, pair_pi_di = self.params["cij"], self.params["tij"], self.params['di'], self.params['pair_pi_di']
        x, v, T_node = self.vars_['x'], self.vars_["v"], self.vars_['T_node']
        z = self.vars_['z']

        if self.variable_substitution: f1 = gb.quicksum(cij[(self.base(i),self.base(j))] * v[i,j] for (i,j) in v.keys())
        else: f1 = gb.quicksum(cij[self.base(i), self.base(j)] * x[k,i,j] for (k,i,j) in x.keys())

        f2 = gb.quicksum(T_node[d] - T_node[p] - di[self.base(p)] - tij[self.base(p),self.base(d)] for p,d in pair_pi_di.items())
        f3 = gb.quicksum((T_node[d] - T_node[p] - di[self.base(p)]) / tij[self.base(p),self.base(d)] for p,d in pair_pi_di.items())

        if self.MoPS:
            D_M = self.sets['D_M']
            if self.variable_substitution:
                f4 = gb.quicksum(v[i,j] for i in D_M for j in N if (i,j) in v.keys())
            else:
                f4 = gb.quicksum(x[k,i,j] for k in K for i in D_M for j in N if (k,i,j) in x.keys())
        else: f4 = 0

        w.append(100)
        f5 = 0

        if self.timetabled_departures:
            f5 = gb.quicksum(
                cij[self.base(i), self.base(j)] * z[d, i, j]
                for i in C
                for j in C
                if (self.base(i), self.base(j)) in self.params['Departures']
                for d in self.params['Departures'][(self.base(i), self.base(j))].keys()
                if (d, i, j) in z  # optional guard, ensures variable exists
            )

        self.m.setObjective(w[0] * f1 + w[1] * f2 + w[2] * f3 - w[3] * f4, gb.GRB.MINIMIZE)

    # === Constraint groups ===
    def add_vehicle_logic_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        v, x, y = self.vars_["v"], self.vars_["x"], self.vars_["y"]
        zeroDepot_node, endDepot_node = zeroDepot, endDepot

        # === CASE 1: VARIABLE SUBSTITUTION ACTIVE ===
        if self.variable_substitution:

            # (1) Each pickup/dropoff visited exactly once
            for i in P + D:
                self.m.addConstr(
                    gb.quicksum(v[i, j] for j in N if (i, j) in v.keys()) == 1,
                    name=f"∑_j v[i,j] == 1  ∀i∈(P∪D)  [outgoing_from_{i}]"
                )
                self.m.addConstr(
                    gb.quicksum(v[j, i] for j in N if (j, i) in v.keys()) == 1,
                    name=f"∑_j v[j,i] == 1  ∀i∈(P∪D)  [incoming_to_{i}]"
                )

            # (2) Each transfer node used at most once
            for i in C:
                self.m.addConstr(
                    gb.quicksum(v[i, j] for j in N if (i, j) in v.keys()) <= 1,
                    name=f"∑_j v[i,j] ≤ 1  ∀i∈C  [transfer_out_{i}]"
                )
                self.m.addConstr(
                    gb.quicksum(v[j, i] for j in N if (j, i) in v.keys()) <= 1,
                    name=f"∑_j v[j,i] ≤ 1  ∀i∈C  [transfer_in_{i}]"
                )

            # (3) Depot departure: each vehicle leaves at most once
            for k in K:
                self.m.addConstr(
                    gb.quicksum(x[k, zeroDepot_node, j] for j in N if (k, zeroDepot_node, j) in x.keys()) <= 1,
                    name=f"∑_j x[k,{zeroDepot_node},j] ≤ 1  ∀k∈K"
                )

            # (4) Link depot departure to vehicle arcs
            for j in N:
                if (zeroDepot_node, j) in v.keys():
                    self.m.addConstr(
                        gb.quicksum(x[k, zeroDepot_node, j] for k in K) == v[zeroDepot_node, j],
                        name=f"∑_k x[k,{zeroDepot_node},j] == v[{zeroDepot_node},{j}]"
                    )

            # (5) Flow conservation at all non-depot nodes
            for i in N:
                if i not in [zeroDepot_node, endDepot_node]:
                    self.m.addConstr(
                        gb.quicksum(v[j, i] for j in N if (j, i) in v.keys())
                        - gb.quicksum(v[i, j] for j in N if (i, j) in v.keys())
                        == 0,
                        name=f"∑_j v[j,i] - ∑_j v[i,j] = 0  ∀i∈(P∪D∪C)"
                    )

            for k in K:
                for i in N:
                    if i not in [zeroDepot_node, endDepot_node]:
                        self.m.addConstr(
                            gb.quicksum(x[k, j, i] for j in N if (k, j, i) in x.keys())
                            - gb.quicksum(x[k, i, j] for j in N if (k, i, j) in x.keys())
                            == 0,
                            name=f"∑_j x[{k},j,i] - ∑_j x[{k},i,j] = 0  ∀i∈(P∪D∪C),∀k"
                        )

            # (6) All vehicles must end at depot
            self.m.addConstr(
                gb.quicksum(v[i, endDepot_node] for i in N if (i, endDepot_node) in v.keys()) == len(K),
                name=f"∑_i v[i,{endDepot_node}] == |K|  [vehicles_end_depot]"
            )

            # (7) Each vehicle ends at depot once
            for k in K:
                self.m.addConstr(
                    gb.quicksum(x[k, i, endDepot_node] for i in N if (k, i, endDepot_node) in x.keys()) == 1,
                    name=f"∑_i x[k,i,{endDepot_node}] == 1  ∀k∈K"
                )

            # (8) Vehicle-depot link consistency
            for i in N:
                if (i, endDepot_node) in v.keys():
                    self.m.addConstr(
                        gb.quicksum(x[k, i, endDepot_node] for k in K) == v[i, endDepot_node],
                        name=f"∑_k x[k,i,{endDepot_node}] == v[i,{endDepot_node}]"
                    )

            # (9) Passenger-vehicle link: y[r,i,j] ≤ v[i,j]
            for r in R:
                for (i, j) in v.keys():
                    if (r, i, j) in y.keys():
                        self.m.addConstr(
                            y[r, i, j] <= v[i, j],
                            name=f"y[{r},{i},{j}] ≤ v[{i},{j}]"
                        )

        # === CASE 2: NO VARIABLE SUBSTITUTION ===
        else:
            # (1') Each pickup/dropoff visited once by some vehicle
            for i in P + D:
                self.m.addConstr(
                    gb.quicksum(x[k, i, j] for k in K for j in N if (k, i, j) in x.keys()) == 1,
                    name=f"∑_k∑_j x[k,{i},j] == 1  ∀i∈(P∪D)"
                )
                # self.m.addConstr(
                #     gb.quicksum(x[k, j, i] for k in K for j in N if (k, j, i) in x.keys()) == 1,
                #     name=f"∑_k∑_j x[k,j,{i}] == 1  ∀i∈(P∪D)"
                # )

            #################################################################################################################################################################

            # A transfer node can be visited at most once
            # for i in C:
            #     self.m.addConstr(
            #         gb.quicksum(x[k, i, j] for k in K for j in N if (i, j) in A) <= 1,
            #         name=f"∑_k∑_j x[k,{i},j] <= 1  ∀i∈(C)"
            #     )
            #     self.m.addConstr(
            #         gb.quicksum(x[k, j, i] for k in K for j in N if (j, i) in A) <= 1,
            #         name=f"∑_k∑_j x[k,j,{i}] <= 1  ∀i∈(C)"
            #     )

            # (2') Each vehicle leaves depot once
            for k in K:
                self.m.addConstr(
                    gb.quicksum(x[k, zeroDepot_node, j] for j in N if (k, zeroDepot_node, j) in x.keys()) == 1,
                    name=f"∑_j x[k,{zeroDepot_node},j] == 1  ∀k∈K"
                )

            # (3') Vehicle flow conservation
            for k in K:
                for i in N:
                    if i not in [zeroDepot_node, endDepot_node]:
                        self.m.addConstr(
                            gb.quicksum(x[k, j, i] for j in N if (k, j, i) in x.keys())
                            - gb.quicksum(x[k, i, j] for j in N if (k, i, j) in x.keys())
                            == 0,
                            name=f"∑_j x[{k},j,i] - ∑_j x[{k},i,j] = 0  ∀i∈(P∪D∪C),∀k"
                        )

            # (4') Each vehicle ends at depot
            for k in K:
                self.m.addConstr(
                    gb.quicksum(x[k, i, endDepot_node] for i in N if (k, i, endDepot_node) in x.keys()) == 1,
                    name=f"∑_i x[k,i,{endDepot_node}] == 1  ∀k∈K"
                )

            # (5') Passenger-vehicle link: y[r,i,j] ≤ ∑_k x[k,i,j]
            for r in R:
                for (i, j) in A:
                    if (r, i, j) in y:
                        self.m.addConstr(
                            y[r, i, j] <= gb.quicksum(x[k, i, j] for k in K),
                            name=f"y[{r},{i},{j}] ≤ ∑_k x[k,{i},{j}]"
                        )

    def add_passenger_balance_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        y, z = self.vars_["y"], self.vars_["z"]
        Cr, fi_r, Departures = self.params["Cr"], self.params["fi_r"], self.params["Departures"]

        if self.timetabled_departures and Departures is not None:

            # (1) Passenger balance for non-transfer nodes
            for r in R:
                for i in P + D:
                    sum_pickup  = gb.quicksum(y[r, i, j] for j in N if (r, i, j) in y.keys())
                    sum_dropoff = gb.quicksum(y[r, j, i] for j in N if (r, j, i) in y.keys())
                    self.m.addConstr(
                        sum_pickup - sum_dropoff == fi_r[r, i],
                        name=f"∑_j y[{r},i,j] - ∑_j y[{r},j,i] = fi_r[{r},{i}]  ∀i∉C"
                    )

            # (2) Passenger balance for transfer nodes with timetabled departures

            # added_cnt = 0
            # empty_rows = []

            for r in R:
                for i in C:
                    sum_pickup  = gb.quicksum(y[r, i, j] for j in N if (r, i, j) in y.keys())
                    sum_dropoff = gb.quicksum(y[r, j, i] for j in N if (r, j, i) in y.keys())

                    z_sum_pickup = gb.quicksum(
                        z[d, r, i, j]
                        for j in C
                        if (self.base(i), self.base(j)) in Departures
                        for d in Departures[(self.base(i), self.base(j))]
                    )


                    z_sum_dropoff = gb.quicksum(
                        z[d, r, j, i]
                        for j in C
                        if (self.base(j), self.base(i)) in Departures
                        for d in Departures[(self.base(j), self.base(i))]
                    )

                    self.m.addConstr(
                        sum_pickup + z_sum_pickup - sum_dropoff - z_sum_dropoff == fi_r[r, i],
                        name=(
                            f"∑_j y[{r},i,j] + ∑_d∑_j z[d,i,j] - ∑_j y[{r},j,i] - ∑_d∑_j z[d,j,i] "
                            f"== fi_r[{r},{i}]  ∀i∈C"
                        )
                    )

                # print(f"[DEBUG] PT passenger-balance constraints added: {added_cnt}")
                # print(f"[DEBUG] Empty rows (would be removed by presolve): {empty_rows[:10]} (total {len(empty_rows)})")

        # === ELSE: Non-timetabled version ===
        else:
            # (1') Passenger balance for non-transfer nodes
            if self.MoPS:
                for r in R:
                    for i in N:
                        if i not in Cr[r]:
                            sum_pickup  = gb.quicksum(y[r, i, j] for j in N if (r, i, j) in y.keys())
                            sum_dropoff = gb.quicksum(y[r, j, i] for j in N if (r, j, i) in y.keys())
                            self.m.addConstr(
                                sum_pickup - sum_dropoff == fi_r[r, i],
                                name=f"∑_j y[{r},i,j] - ∑_j y[{r},j,i] = fi_r[{r},{i}]  ∀i∉C"
                            )
            
            else:
                for r in R:
                    for i in N:
                        if i not in Cr[r]:
                            sum_pickup  = gb.quicksum(y[r, i, j] for j in N if (r, i, j) in y.keys())
                            sum_dropoff = gb.quicksum(y[r, j, i] for j in N if (r, j, i) in y.keys())
                            self.m.addConstr(
                                sum_pickup - sum_dropoff == fi_r[r, i],
                                name=f"∑_j y[{r},i,j] - ∑_j y[{r},j,i] = fi_r[{r},{i}]  ∀i∉C"
                            )

            # (2') Passenger balance for transfer nodes (without timetabled departures)
            if Cr:
                for r in R:
                    for i in Cr[r]:
                        sum_pickup  = gb.quicksum(y[r, i, j] for j in N if (r, i, j) in y.keys())
                        sum_dropoff = gb.quicksum(y[r, j, i] for j in N if (r, j, i) in y.keys())

                        z_sum_pickup = gb.quicksum(
                            z[(i, j)] for j in Cr[r] if (i, j) in z.keys()
                        )
                        z_sum_dropoff = gb.quicksum(
                            z[(j, i)] for j in Cr[r] if (j, i) in z.keys()
                        )

                        self.m.addConstr(
                            sum_pickup + z_sum_pickup - sum_dropoff - z_sum_dropoff == fi_r[r, i],
                            name=(
                                f"∑_j y[{r},i,j] + ∑_j z[i,j] - ∑_j y[{r},j,i] - ∑_j z[j,i] "
                                f"= fi_r[{r},{i}]  ∀i∈Cr[{r}]"
                            )
                        )

    def add_capacity_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, y = self.vars_['x'], self.vars_["v"], self.vars_["y"]
        qr, Q = self.params["qr"], self.params["Q"]
        # === CASE 1: Variable Substitution Active ===
        if self.variable_substitution:
            for (i, j) in v.keys():
                self.m.addConstr(
                    gb.quicksum(qr[r] * y[r, i, j] for r in R if (r, i, j) in y.keys()) <= Q * v[i, j],
                    name=(
                        f"∑_r (q_r[r]*y[r,{i},{j}]) ≤ Q * v[{i},{j}]  ∀(i,j)∈A"
                    )
                )

        # === CASE 2: Classic x[k,i,j]-based model ===
        else:
            for (i, j) in A:
                if not (i in C and j in C):
                    self.m.addConstr(
                        gb.quicksum(qr[r] * y[r, i, j] for r in R if (r, i, j) in y.keys())
                        <= Q * gb.quicksum(x[k, i, j] for k in K),
                        name=(
                            f"∑_r (q_r[r]*y[r,{i},{j}]) ≤ Q * ∑_k x[k,{i},{j}]  ∀(i,j)∈A"
                        )
                    )

    def add_time_consistency_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, z, T_node, T_veh, D_node = self.vars_["x"], self.vars_["v"], self.vars_["z"], self.vars_["T_node"], self.vars_["T_veh"], self.vars_["D_node"]
        tij, di, ei, li, ek, lk, Lbar, M, T, Cr = self.params["tij"], self.params["di"], self.params["ei"], self.params["li"], self.params["ek"], self.params["lk"], self.params['Lbar'], self.params["M"], self.params["T"], self.params['Cr']
        pair_pi_di = self.params["pair_pi_di"]
        zeroDepot_node, endDepot_node = zeroDepot, endDepot

        if self.variable_substitution:
            # (1) Depot → First node timing
            for k in K:
                for j in N:
                    if (k, zeroDepot_node, j) in x.keys():
                        self.m.addConstr(
                            T_node[j] >= T_veh[k] + tij[(self.base(zeroDepot_node), self.base(j))] - M*(1 - x[k, zeroDepot_node, j]),
                            name=f"T[{j}] ≥ Tveh[{k}] + t[{self.base(zeroDepot_node)},{self.base(j)}] - M(1 - x[{k},{zeroDepot_node},{j}])"
                        )

            # (2) Service time progression between served nodes
            for (i, j) in v.keys():
                self.m.addConstr(
                    T_node[j] >= T_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(j))] - M*(1 - v[i, j]),
                    name=f"T[{j}] ≥ T[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(j)}] - M(1 - v[{i},{j}])"
                )

            # (3) Ride-time constraints between pickup–dropoff pairs
            for p, d in pair_pi_di.items():
                self.m.addConstr(
                    T_node[d] - (T_node[p] + di[self.base(p)]) >= tij[(self.base(p), self.base(d))],
                    name=f"T[{d}] - (T[{p}] + d[{self.base(p)}]) ≥ t[{self.base(p)},{self.base(d)}]"
                )
                self.m.addConstr(
                    T_node[d] - (T_node[p] + di[self.base(p)]) <= Lbar[self.base(p)],
                    name=f"T[{d}] - (T[{p}] + d[{self.base(p)}]) ≤ L̄[{self.base(p)}]"
                )

            # (4) Time windows at pickup/dropoff nodes
            for i in P + D:
                self.m.addConstr(
                    T_node[i] >= ei[self.base(i)],
                    name=f"T[{i}] ≥ e[{self.base(i)}]"
                )
                self.m.addConstr(
                    T_node[i] <= li[self.base(i)],
                    name=f"T[{i}] ≤ l[{self.base(i)}]"
                )

            # (5) Cumulative time between service nodes (duration tracking)
            for (i, j) in v.keys():
                self.m.addConstr(
                    D_node[j] >= D_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(j))] - M*(1 - v[i, j]),
                    name=f"D[{j}] ≥ D[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(j)}] - M(1 - v[{i},{j}])"
                )

            # (6) Maximum route duration
            for i in N:
                if (i, endDepot_node) in A:
                    self.m.addConstr(
                        D_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(endDepot_node))] <= T,
                        name=f"D[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(endDepot_node)}] ≤ T"
                    )

            # (7) Depot departure window
            for k in K:
                self.m.addConstr(
                    T_veh[k] >= ek[k],
                    name=f"Tveh[{k}] ≥ e_k[{k}]"
                )
                self.m.addConstr(
                    T_veh[k] <= lk[k],
                    name=f"Tveh[{k}] ≤ l_k[{k}]"
                )

            # (8) Transfer time links (fixed-route for Cr)
            if Cr:
                for r in R:
                    for i in Cr[r]:
                        for j in Cr[r]:
                            if (self.base(i), self.base(j)) in tij:
                                self.m.addConstr(
                                    T_node[j] >= T_node[i] + tij[(self.base(i), self.base(j))] - M*(1 - z[(i, j)]),
                                    name=f"T[{j}] ≥ T[{i}] + t[{self.base(i)},{self.base(j)}] - M(1 - z[{i},{j}])"
                                )

        else:
            # (9) Depot → First node (classic x formulation)
            for k in K:
                for j in N:
                    if (k, zeroDepot_node, j) in x.keys():
                        self.m.addConstr(
                            T_node[j] >= T_veh[k] + tij[(self.base(zeroDepot_node), self.base(j))] - M*(1 - x[k, zeroDepot_node, j]),
                            name=f"T[{j}] ≥ Tveh[{k}] + t[{self.base(zeroDepot_node)},{self.base(j)}] - M(1 - x[{k},{zeroDepot_node},{j}])"
                        )

            # (10) Time consistency between service nodes
            for (i, j) in A:
                # if not (i in C and j in C):
                self.m.addConstr(
                    T_node[j] >= T_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(j))] - M*(1 - gb.quicksum(x[k, i, j] for k in K if (k, i, j) in x.keys())),
                    name=f"T[{j}] ≥ T[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(j)}] - M(1 - ∑_k x[k,{i},{j}])"
                )

            ######################################################################################################################################################

            # (11) Transfer route timing for Cr requests
            for r in R:
                if Cr:
                    for i in Cr[r]:
                        for j in Cr[r]:
                            if (i, j) in A:
                                self.m.addConstr(
                                    T_node[j] >= T_node[i] + tij[(self.base(i), self.base(j))] - M*(1 - z[(i, j)]),
                                    name=f"T[{j}] ≥ T[{i}] + t[{self.base(i)},{self.base(j)}] - M(1 - z[{i},{j}])"
                                )

            # (12) Time windows
            for i in P + D:
                self.m.addConstr(
                    T_node[i] >= ei[self.base(i)],
                    name=f"T[{i}] ≥ e[{self.base(i)}]"
                )
                self.m.addConstr(
                    T_node[i] <= li[self.base(i)],
                    name=f"T[{i}] ≤ l[{self.base(i)}]"
                )

            # (13) Ride-time constraints
            for i in P:
                drop = pair_pi_di[i]
                if (i, drop) in A:
                    self.m.addConstr(
                        T_node[drop] - (T_node[i] + di[self.base(i)]) >= tij[(self.base(i), self.base(drop))],
                        name=f"T[{drop}] - (T[{i}] + d[{self.base(i)}]) ≥ t[{self.base(i)},{self.base(drop)}]"
                    )
                    self.m.addConstr(
                        T_node[drop] - (T_node[i] + di[self.base(i)]) <= Lbar[self.base(i)],
                        name=f"T[{drop}] - (T[{i}] + d[{self.base(i)}]) ≤ L̄[{self.base(i)}]"
                    )

            # (14) Cumulative duration D tracking
            # for (i, j) in A:
            #     if not (i in C and j in C):
            #         self.m.addConstr(
            #             D_node[j] >= D_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(j))] - M*(1 - gb.quicksum(x[k, i, j] for k in K)),
            #             name=f"D[{j}] ≥ D[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(j)}] - M(1 - ∑_k x[k,{i},{j}])"
            #         )

            for i in N:
                for j in N:
                    if (i,j) in A:
                        self.m.addConstr(
                            D_node[j] >= D_node[i] + di[self.base(i)] + tij[self.base(i),self.base(j)] - M*(1 - gb.quicksum(x[k,i,j] for k in K if (k, i, j) in x.keys())),
                            name=f"D_link[{i},{j}]"
                        )

            #################################################################################################################################################################################################

            # (15) Route duration limit
            for i in N:
                if (i, endDepot_node) in A:
                    self.m.addConstr(
                        D_node[i] + di[self.base(i)] + tij[(self.base(i), self.base(endDepot_node))] <= T,
                        name=f"D[{i}] + d[{self.base(i)}] + t[{self.base(i)},{self.base(endDepot_node)}] ≤ T"
                    )

            # (16) Vehicle departure windows
            for k in K:
                self.m.addConstr(
                    T_veh[k] >= ek[k],
                    name=f"Tveh[{k}] ≥ e_k[{k}]"
                )
                self.m.addConstr(
                    T_veh[k] <= lk[k],
                    name=f"Tveh[{k}] ≤ l_k[{k}]"
                )

    def add_variable_substitution_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        v, x = self.vars_["v"], self.vars_["x"]
        zeroDepot_node = zeroDepot

        # for (i, j) in A:
        #     if not (i in C and j in C):
        #         self.m.addConstr(
        #             v[i, j] == gb.quicksum(x[k, i, j] for k in K),
        #             name=f"v[{i},{j}] = ∑_k x[k,{i},{j}]  ∀(i,j)∈A, i≠{zeroDepot_node}"
        #         )

        for (i,j) in A:
            if i != zeroDepot_node:
                if (i,j) in v.keys():
                    self.m.addConstr(v[i,j] == gb.quicksum(x[k,i,j] for k in K))

        #####################################################################################################################################################################################################

    def add_transfer_node_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, z = self.vars_['x'], self.vars_["v"], self.vars_["z"]
        Cr, Departures = self.params['Cr'], self.params["Departures"]

        if self.variable_substitution:
            # === CASE: Variable Substitution with Timetabled Departures ===
            if self.timetabled_departures and Departures is not None:
                skip = 1
                # for i in C:
                #     self.m.addConstr(
                #         gb.quicksum(v[j, i] for j in N if (j, i) in A and not j in C)
                #         ==
                #         gb.quicksum(
                #             z[(d, i, j)] + z[(d, j, i)]
                #             for j in C
                #             if (i, j) in A and i != j and (self.base(i), self.base(j)) in Departures
                #             for d in Departures[(self.base(i), self.base(j))].keys()
                #         ),
                #         name=(
                #             f"∑_j v[j,{i}] = ∑_d∑_j (z[d,{i},j] + z[d,j,{i}]) "
                #             f"∀i∈C with timetabled Departures"
                #         )
                #     )

            # === CASE: Variable Substitution without Timetabled Departures ===
            # else:
            #     if Cr:
            #         for r in R:
            #             for i in Cr[r]:
            #                 self.m.addConstr(
            #                     gb.quicksum(v[j, i] for j in N if (j, i) in A and not j in C)
            #                     ==
            #                     gb.quicksum(
            #                         z[(i, j)] + z[(j, i)]
            #                         for j in Cr[r] if (i, j) in A
            #                     ),
            #                     name=(
            #                         f"∑_j v[j,{i}] = ∑_j (z[{i},j] + z[j,{i}]) "
            #                         f"∀i∈Cr[{r}], variable_substitution"
            #                     )
            #                 )

            #         for r in R:
            #             self.m.addConstr(
            #                 gb.quicksum(v[i, j] for i in Cr[r] for j in N if (i, j) in A and not (i in C and j in C)) <= 2,
            #                 name=f"∑_i∑_j v[i,j] ≤ 2  ∀i∈Cr[{r}]"
            #             )

            else:
                if Cr:
                    for r in R:
                        for i in Cr[r]:
                            self.m.addConstr(
                                gb.quicksum(v[j,i] for j in N if (j,i) in v.keys()) ==
                                gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) and (j,i) in z.keys())
                            )
                    for r in R:
                        self.m.addConstr(gb.quicksum(v[i,j] for i in Cr[r] for j in N if (i,j) in v.keys()) <= 2)

        # === CASE: No Variable Substitution ===
        else:
            if self.timetabled_departures and Departures is not None:
                skip = 1
                # for i in C:
                #     self.m.addConstr(
                #         gb.quicksum(x[k, j, i] for k in K for j in N if (k, j, i) in x.keys())
                #         ==
                #         gb.quicksum(
                #             z[(d, i, j)] + z[(d, j, i)]
                #             for j in C
                #             if (i, j) in A and i != j and (self.base(i), self.base(j)) in Departures
                #             for d in Departures[(self.base(i), self.base(j))].keys()
                #         ),
                #         name=(
                #             f"∑_j v[j,{i}] = ∑_d∑_j (z[d,{i},j] + z[d,j,{i}]) "
                #             f"∀i∈C with timetabled Departures"
                #         )
                #     )

            else:
                if Cr:
                    for r in R:
                        for i in Cr[r]:
                            self.m.addConstr(
                                gb.quicksum(x[k, j, i] for k in K
                                    for j in N if (k, j, i) in x.keys()
                                )
                                ==
                                gb.quicksum(
                                    z[(i, j)] + z[(j, i)]
                                    for j in Cr[r] if (i, j) in A
                                ),
                                name=(
                                    f"∑_k∑_j x[k,j,{i}] = ∑_j (z[{i},j] + z[j,{i}]) "
                                    f"∀i∈Cr[{r}], no_variable_substitution"
                                )
                            )

    def add_battery_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, B_node, E_node = self.vars_['x'], self.vars_["v"], self.vars_["B_node"], self.vars_["E_node"]
        tij = self.params["tij"]
        beta, alpha, C_bat, gamma = self.params["beta"], self.params["alpha"], self.params["C_bat"], self.params["gamma"]

        if self.variable_substitution:
            # (1) Battery SOC variation along arcs (with variable substitution)
            for (i, j) in v.keys():
                if i not in F and j not in F:
                    self.m.addConstr(
                        B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] - C_bat * (1 - v[i, j]) <= 0,
                        name=f"B[{j}] - B[{i}] + β * t[{self.base(i)},{self.base(j)}] - Cbat * (1 - v[{i},{j}]) ≤ 0"
                    )
                    self.m.addConstr(
                        B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] + C_bat * (1 - v[i, j]) >= 0,
                        name=f"B[{j}] - B[{i}] + β * t[{self.base(i)},{self.base(j)}] + Cbat * (1 - v[{i},{j}]) ≥ 0"
                    )

        else:
            # (1') Battery SOC variation along arcs (without variable substitution)
            for (i, j) in A:
                if i not in F and j not in F:
                    self.m.addConstr(
                        B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] - C_bat * (1 - gb.quicksum(x[k, i, j] for k in K if (k, i, j) in x.keys())) <= 0,
                        name=f"B[{j}] - B[{i}] + β * t[{self.base(i)},{self.base(j)}] - Cbat * (1 - ∑_k x[k,{i},{j}]) ≤ 0"
                    )

        # (2) Battery capacity constraint at charging stations
        for s in F:
            self.m.addConstr(
                C_bat - (B_node[s] + alpha * E_node[s]) >= 0,
                name=f"Cbat - (B[{s}] + α * E[{s}]) ≥ 0  ∀s∈F"
            )

        # (3) Minimum state of charge constraint
        for i in N:
            self.m.addConstr(
                B_node[i] >= gamma * C_bat,
                name=f"B[{i}] ≥ γ * Cbat  ∀i∈N"
            )

        # (4) Initial state of charge at depot
        self.m.addConstr(
            B_node[zeroDepot] == C_bat,
            name=f"B[{zeroDepot}] = Cbat  [Initial SOC]"
)

    def add_scheduled_PT_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        z, a, T_node, y = self.vars_["z"], self.vars_["a"], self.vars_["T_node"], self.vars_["y"]
        M, tij, Departures = self.params["M"], self.params['tij'], self.params["Departures"]

        # (1) PT departure occurs after service start at transfer node i
        for r in R:
            for i in C:
                expr = (
                    M + gb.quicksum(
                        (Departures[self.base(i), self.base(j)][d] - M) * z[d, r, i, j]
                        for j in C
                        if (i,j) in A
                        for d in Departures[self.base(i), self.base(j)]
                    ) - T_node[i]
                )
                self.m.addConstr(
                    expr >= 0,
                    name=f"M + ∑_d∑_j (Dep[{self.base(i)},j][d] - M)·z[d,{i},j] - T[{i}] ≥ 0  [PT_depart_after_service]"
                )

        # (2) Service at arrival transfer node j starts after PT arrival
        for r in R:
            for j in C:
                expr = (
                    T_node[j]
                    - gb.quicksum(
                        (Departures[self.base(i), self.base(j)][d] + tij[self.base(i), self.base(j)]) * z[d, r, i, j]
                        for i in C
                        if (i,j) in A
                        for d in Departures[self.base(i), self.base(j)]
                    )
                )
                self.m.addConstr(
                    expr >= 0,
                    name=f"T[{j}] - ∑_d∑_i (Dep[{self.base(i)},j][d] + t[{self.base(i)},j])·z[d,{i},j] ≥ 0  [PT_arrival_before_service]"
                )

        # (3) Node visit upper bound: if used in PT trip → mark visited
        for i in C:
            expr_lhs = M * a[i]
            expr_rhs = (
                gb.quicksum(
                    z[d, r, i, j]
                    for j in C
                    if (i,j) in A
                    for d in Departures[self.base(i), self.base(j)]
                    for r in R
                )
                + gb.quicksum(
                    z[d, r, j, i]
                    for j in C
                    if (i,j) in A
                    for d in Departures[self.base(j), self.base(i)]
                    for r in R
                )
            )
            self.m.addConstr(
                expr_lhs >= expr_rhs,
                name=f"M·a[{i}] ≥ ∑_d∑_j∑_r z[d,r,{i},j] + ∑_d∑_j∑_r z[d,r,j,{i}]  [PT_node_visit_upper]"
            )

        # (4) Node visit lower bound: node marked as visited only if used
        for i in C:
            expr_rhs = (
                gb.quicksum(
                    z[d, r, i, j]
                    for j in C
                    if (i,j) in A
                    for d in Departures[self.base(i), self.base(j)]
                    for r in R
                )
                + gb.quicksum(
                    z[d, r, j, i]
                    for j in C
                    if (i,j) in A
                    for d in Departures[self.base(j), self.base(i)]
                    for r in R
                )
            )
            self.m.addConstr(
                a[i] <= expr_rhs,
                name=f"a[{i}] ≤ ∑_d∑_j z[d,{i},j] + ∑_d∑_j z[d,j,{i}]  [PT_node_visit_lower]"
            )

    def add_artificial_node_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        M = self.params['M']
        x, v, z, a, T_node = self.vars_['x'], self.vars_['v'], self.vars_["z"], self.vars_["a"], self.vars_['T_node']
        for (i, m_) in N:
            if m_ < max(n for (j, n) in N if j == i):
                # (1) Monotonic activation of layers/modes
                self.m.addConstr(
                    a[i, m_] - a[i, m_ + 1] >= 0,
                    name=f"a[{i},{m_}] - a[{i},{m_+1}] ≥ 0  [Monotonic_activation_i={i}]"
                )

        # for j in N:
        #     if self.variable_substitution:
        #         self.m.addConstr(a[j] - gb.quicksum(v[i,j] for i in N if (i,j) in A) == 0 , name = "link a_im to v_ij")

        #     else:
        #         self.m.addConstr(a[j] - gb.quicksum(x[k, i, j] for k in K for i in N if (i,j) in A) == 0 , name = "link a_im to v_ij")

    def add_MoPS_constraints(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v = self.vars_['x'], self.vars_['v']
        pair_pi_di_M = self.params['pair_pi_di_M']
        # (1) MoPS entry and exit limits
        if self.variable_substitution:
            for i in P_M + D_M:
                # Exit constraint
                self.m.addConstr(
                    gb.quicksum(v[i, j] for j in N if (i, j) in A) <= 1,
                    name=f"∑_j v[{i},j] ≤ 1  ∀i∈(P_M∪D_M)  [MoPS_exit_limit_i={i}]"
                )

                # Entry constraint
                self.m.addConstr(
                    gb.quicksum(v[j, i] for j in N if (j, i) in A) <= 1,
                    name=f"∑_j v[j,{i}] ≤ 1  ∀i∈(P_M∪D_M)  [MoPS_entry_limit_i={i}]"
                )

            # (2) Each picked-up MoPS request must be dropped off
            for i in P_M:
                self.m.addConstr(
                    gb.quicksum(v[i, j] for j in N if (i, j) in A)
                    - gb.quicksum(v[pair_pi_di_M[i], j] for j in N if (pair_pi_di_M[i], j) in A)
                    == 0,
                    name=(
                        f"∑_j v[{i},j] - ∑_j v[{pair_pi_di_M[i]},j] = 0  "
                        f"∀i∈P_M  [MoPS_pickup_drop_link_i={i}]"
                    )
                )

        else:
            for i in P_M + D_M:
                # Exit constraint
                self.m.addConstr(
                    gb.quicksum(x[k, i, j] for k in K for j in N if (k, i, j) in x.keys()) <= 1,
                    name=f"∑_j v[{i},j] ≤ 1  ∀i∈(P_M∪D_M)  [MoPS_exit_limit_i={i}]"
                )

                # Entry constraint
                self.m.addConstr(
                    gb.quicksum(x[k, j, i] for k in K for j in N if (k, j, i) in x.keys()) <= 1,
                    name=f"∑_j v[j,{i}] ≤ 1  ∀i∈(P_M∪D_M)  [MoPS_entry_limit_i={i}]"
                )

            # (2) Each picked-up MoPS request must be dropped off
            for i in P_M:
                self.m.addConstr(
                    gb.quicksum(x[k, i, j] for k in K for j in N if (k, i, j) in x.keys())
                    - gb.quicksum(x[k, pair_pi_di_M[i], j] for k in K for j in N if (k, pair_pi_di_M[i], j) in x.keys())
                    == 0,
                    name=(
                        f"∑_j v[{i},j] - ∑_j v[{pair_pi_di_M[i]},j] = 0  "
                        f"∀i∈P_M  [MoPS_pickup_drop_link_i={i}]"
                    )
                )

    def add_imjn_enhancements(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets['P_M'], self.sets['D_M'], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        v, x, y = self.vars_["v"], self.vars_["x"], self.vars_["y"]
        zeroDepot_node, endDepot_node = zeroDepot, endDepot

        for (i,j) in x:
            if (j,i) in A:
                if self.variable_substitution:
                    self.m.addConstr(v[i,j] + v[j,i] <= 1, name = "A vehicle can only travel in one direction")
                
                else:
                    self.m.addConstr(gb.quicksum(x[k,i,j] + x[k,j,i] for k in K) <= 1, name = "A vehicle can only travel in one direction")

    def base_model_optimal_solution_constrainer_paper(self):
        nodes, N, P, D, C, F, R, K, P_M, D_M, zeroDepot, endDepot, A = (
            self.sets["nodes"],
            self.sets["N"],
            self.sets["P"],
            self.sets["D"],
            self.sets["C"],
            self.sets["F"],
            self.sets["R"],
            self.sets["K"],
            self.sets['P_M'],
            self.sets['D_M'],
            self.sets["zeroDepot"],
            self.sets["endDepot"],
            self.sets["A"],
        )
        v, x, y, z = self.vars_["v"], self.vars_["x"], self.vars_["y"], self.vars_["z"]

        if self.use_imjn:
            # -------------------------------
            # IMJN version (layered nodes)
            # -------------------------------

            # Vehicle 1
            # vehicle1_arcs = []
            vehicle1_arcs = [
                ((0, 1), (1, 1)),
                ((1, 1), (3, 1)),
                ((3, 1), (2, 1)),
                ((2, 1), (5, 1)),
                ((5, 1), (10,1)),
                ((10,1), (10,2)),
                ((10,2), (8, 1)),
                ((8, 1), (9, 1)),
            ]

            # Vehicle 2
            vehicle2_arcs = [
                ((0, 1), (11,1)),
                ((11,1), (7, 1)),
                ((7, 1), (4, 1)),
                ((4, 1), (12,1)),
                ((12,1), (6, 1)),
                ((6, 1), (9, 1)),
            ]

            # # Requests
            # request_arcs = []
            request_arcs = {
                1: [((1, 1), (3, 1)), ((3, 1), (2, 1)), ((2, 1), (5, 1))],
                2: [((2, 1), (5, 1)), ((5, 1), (10, 1)), ((11, 1), (7, 1)), ((7, 1), (4, 1)), ((4, 1), (12, 1)), ((12, 1), (6, 1))],
                3: [((3, 1), (2, 1)), ((2, 1), (5, 1)), ((5, 1), (10, 1)), ((11, 1), (7, 1))],
                4: [((4, 1), (12, 1)), ((10, 2), (8, 1))],
            }

            transfer_arcs = [
                ((10, 1), (11, 1)),
                ((12, 1), (10, 2)),
            ]
            # transfer_arcs = []

        else:
            # -------------------------------
            # Non-IMJN version (simple nodes)
            # -------------------------------
            if self.MoPS:
                vehicle1_arcs = [
                    (0, 1),
                    (1, 3),
                    (3, 2),
                    (2, 5),
                    (5, 15),
                    (15, 18),
                    (18, 21),
                    (21, 8),
                    (8, 11),
                ]

                # Vehicle 2
                vehicle2_arcs = [
                    (0, 9),
                    (9, 10),
                    (10, 16),
                    (16, 19),
                    (19, 7),
                    (7, 4),
                    (4, 6),
                    (6, 22),
                    (22, 11),
                ]

                # Requests
                request_arcs = {
                    1: [(1, 3), (3, 2), (2, 5)],
                    2: [(2, 5), (5, 15), (16, 19), (19, 7), (7, 4), (4, 6)],
                    3: [(3, 2), (2, 5), (5, 15), (15, 18), (19, 7)],
                    4: [(4, 6), (6, 22), (21, 8)],
                    5: [(9,10)]
                }

                transfer_arcs = [
                    (15, 16),
                    (18, 19),
                    (22, 21)
                ]

            else:

            # Vehicle 1
                vehicle1_arcs = [
                    (0, 1),

                    (1, 3),
                    (3, 2),
                    (2, 5),
                    (5, 13),
                    (13, 16),
                    (16, 19),
                    (19, 8),
                    (8, 9),
                ]

                # Vehicle 2
                vehicle2_arcs = [
                    (0, 14),
                    (14, 17),
                    (17, 7),
                    (7, 4),
                    (4, 6),
                    (6, 20),
                    (20, 9),
                ]

                # Requests
                request_arcs = {
                    1: [(1, 3), (3, 2), (2, 5)],
                    2: [(2, 5), (5, 13), (14, 17), (17, 7), (7, 4), (4, 6)],
                    3: [(3, 2), (2, 5), (5, 13), (13, 16), (17, 7)],
                    4: [(4, 6), (6, 20), (19, 8)],
                }

                transfer_arcs = [
                    (13, 14),
                    (16, 17),
                    (20, 19)
                ]

        # --- Fix x ---
        for k, arcs in enumerate([vehicle1_arcs, vehicle2_arcs]):
            for (i, j) in arcs:
                self.m.addConstr(x[k, i, j] == 1, name=f"fix_x[{k},{i},{j}]")

        # for k in [0, 1]:
        #     for (i, j) in A:
        #         if (i, j) not in vehicle1_arcs + vehicle2_arcs:
        #             if (k, i, j) in x:
        #                 self.m.addConstr(x[k, i, j] == 0, name=f"fix_x0[{k},{i},{j}]")

        # --- Fix y ---
        for r, arcs in request_arcs.items():
            for (i, j) in arcs:
                if (r, i, j) in y:
                    self.m.addConstr(y[r, i, j] == 1, name=f"fix_y[{r},{i},{j}]")

        for r in R:
            for (i, j) in A:
                if (i, j) not in [a for arcs in request_arcs.values() for a in arcs]:
                    if (r, i, j) in y:
                        self.m.addConstr(y[r, i, j] == 0, name=f"fix_y0[{r},{i},{j}]")


        # if self.timetabled_departures:
        #     temp = 1
        #     for (i, j) in transfer_arcs:
        #         if (i, j) in z:
        #             self.m.addConstr(gb.quicksum(z[d, i, j] for d in self.params['Departures'][(self.base(i), self.base(j))]) >= 1, name=f"fix_z[{i},{j}]")

        #     for (i, j) in A:
        #         if (i, j) not in transfer_arcs:
        #             for d in self.params['Departures'][(self.base(i), self.base(j))]:
        #                 self.m.addConstr(z[d, i, j] == 0, name=f"fix_z0[{i},{j}]")

        # else:
        #     # --- Fix z ---
        #     for (i, j) in transfer_arcs:
        #         if (i, j) in z:
        #             self.m.addConstr(z[i, j] == 1, name=f"fix_z[{i},{j}]")

        #     for (i, j) in A:
        #         if (i, j) not in transfer_arcs:
        #             if (i, j) in z:
        #                 self.m.addConstr(z[i, j] == 0, name=f"fix_z0[{i},{j}]")

    def debugging_constraint(self):
        z = self.vars_['z']
        self.m.addConstr(gb.quicksum(z[d, i, j] for (d, i, j) in z.keys()) >= 1)
