##### Imports ##### 
import gurobipy as gb

class DARPConstraintBuilder:
    def __init__(self, m, vars_, sets, params, duplicate_transfers=True, arc_elimination=True, variable_substitution=True, subtour_elimination=True, transfer_node_strengthening=True,
                 ev_constraints=False, timetabled_departures=False, use_imjn=False, MoPS = False):
        self.duplicate_transfers = duplicate_transfers,
        self.arc_elimination = arc_elimination,
        self.variable_substitution = variable_substitution,
        self.subtour_elimination = subtour_elimination,
        self.transfer_node_strengthening = transfer_node_strengthening,
        self.ev_constraints = ev_constraints,
        self.timetabled_departures = timetabled_departures,
        self.use_imjn = use_imjn
        self.MoPS = MoPS
        self.m = m
        self.vers_= vars_
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

        Departures = self.params.get("Departures", None)

        ### --- Original variables ---
        # Vehicle routing arc variables
        x = self.m.addVars(K, A, vtype=gb.GRB.BINARY, name="x")   # vehicle arcs
        # Passenger flow
        y = self.m.addVars(R, A, vtype=gb.GRB.CONTINUOUS, lb=0, name="y")
        # Node timing variables
        T_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")
        T_veh  = self.m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")
        D_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

        ### --- Developments ---
        # Optional variables
        v, a, B_node, E_node = None, None, None, None
        z = {}

        # Variable Substitution
        if self.variable_substitution:
            v = self.m.addVars(A, vtype=gb.GRB.BINARY, name="v")      # DAR arcs

        # Transfer / fixed-route arcs
        z = {}
        if self.timetabled_departures and Departures is not None:
            # Timetabled departures case
            for (i, j) in A:
                if i in C and j in C and (self.base(i), self.base(j)) in Departures:
                    for d in Departures[(self.base(i), self.base(j))]:
                        z[(d, i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{d},{i},{j}]")
                # Artificial node marker (used in some strengthened formulations)

        if self.use_imjn:
            a = self.m.addVars(N, vtype=gb.GRB.BINARY, name="a")

        elif C:  
            # Transfer Original Modeling
            for i in C:
                for j in C:
                    z[(i, j)] = self.m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")

        # EV-specific variables
        if self.ev_constraints:
            B_node = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")  # battery SOC
            E_node = self.m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")  # recharge time

        return {
            "x": x, "v": v, "y": y, "z": z,
            "T_node": T_node, "T_veh": T_veh, "D_node": D_node,
            "a": a, "B_node": B_node, "E_node": E_node
        }

    # === Objective ===
    def set_objective(self, w = [1,0,0,0]):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        cij, tij, di, pair_pi_di = self.params["cij"], self.params["tij"], self.params['di'], self.params['pair_pi_di']
        x, v, T_node = self.vars_['x'], self.vars_["v"], self.vars_['T_node']
        
        f2 = gb.quicksum(T_node[d] - T_node[p] - di[self.base(p)] - tij[self.base(p),self.base(d)] for p,d in pair_pi_di.items())
        f3 = gb.quicksum((T_node[d] - T_node[p] - di[self.base(p)]) / tij[self.base(p),self.base(d)] for p,d in pair_pi_di.items())

        if self.variable_substitution: f1 = gb.quicksum(cij[(self.base(i),self.base(j))] * v[i,j] for (i,j) in A)

        else: f1 = gb.quicksum(cij[self.base(i), self.base(j)] * x[k,i,j] for k in K for (i,j) in A)

        if self.MoPS:
            D_M = self.sets['D_M']
            f4 = gb.quicksum(v[i,j] for i in D_M for j in N if (i,j) in A)
        else: f4 = 0

        self.m.setObjective(w[0] * f1 + w[1] * f2 + w[2] * f3 - w[3] * f4 , gb.GRB.MINIMIZE)

    # === Constraint groups ===
    def add_vehicle_logic_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        v, x, y = self.vars_["v"], self.vars_["x"], self.vars_["y"]
        zeroDepot_node, endDepot_node = zeroDepot, endDepot

        if self.variable_substitution:
            # pickups/dropoffs visited exactly once
            for i in P + D:
                self.m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in A) == 1)
                self.m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) == 1)

            # transfer nodes at most once
            for i in C:
                self.m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in A) <= 1)
                self.m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) <= 1)

            # depot logic
            for k in K:
                self.m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for j in N if (zeroDepot_node,j) in A) <= 1)
            for j in N:
                if (zeroDepot_node, j) in A:
                    self.m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for k in K) == v[zeroDepot_node,j])

            # flow conservation
            for i in (P+D+C):
                self.m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) -
                            gb.quicksum(v[i,j] for j in N if (i,j) in A) == 0)

            # vehicle must end at depot
            self.m.addConstr(gb.quicksum(v[i,endDepot_node] for i in N if (i,endDepot_node) in A) == len(K))
            for k in K:
                self.m.addConstr(gb.quicksum(x[k,i,endDepot_node] for i in N if (i,endDepot_node) in A) == 1)

            for i in N:
                if (i,endDepot_node) in A:
                    self.m.addConstr(
                        gb.quicksum(x[k,i,endDepot_node] for k in K) == v[i,endDepot_node],
                        name=f"endDepot_v_link[{i}]"
                    )

            # added by me
            for r in R:
                for (i, j) in A:
                    if (r, i, j) in y:
                        # y[r,i,j] ≤ v[i,j]
                        self.m.addConstr(y[r, i, j] <= v[i, j],
                                    name=f"passenger_vehicle_link[{r},{i},{j}]")

        else: 
            # pickups/dropoffs visited exactly once by some vehicle
            for i in P + D:
                self.m.addConstr(gb.quicksum(x[k,i,j] for k in K for j in N if (i,j) in A) == 1)

            # each vehicle leaves depot exactly once
            for k in K:
                self.m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for j in N if (zeroDepot_node,j) in A) == 1)

            # vehicle flow conservation
            for k in K:
                for i in (P + D + C):
                    self.m.addConstr(gb.quicksum(x[k,j,i] for j in N if (j,i) in A) -
                                gb.quicksum(x[k,i,j] for j in N if (i,j) in A) == 0)

            # each vehicle ends at depot
            for k in K:
                self.m.addConstr(gb.quicksum(x[k,i,endDepot_node] for i in N if (i,endDepot_node) in A) == 1)

            # added by me
            for r in R:
                for (i, j) in A:
                    if (r, i, j) in y:
                        # y[r,i,j] ≤ v[i,j]
                        self.m.addConstr(y[r, i, j] <= gb.quicksum(x[k, i, j] for k in K),
                                    name=f"passenger_vehicle_link[{r},{i},{j}]")

    def add_passenger_balance_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        y, z = self.vars_["y"], self.vars_["z"]
        Cr, fi_r, Departures = self.params["Cr"], self.params["fi_r"], self.params["Departures"]

        if self.timetabled_departures and Departures != None:
            for r in R:
                for i in C:
                    lhs_out = gb.quicksum(y[r, i, j] for j in N if (i, j) in A)
                    lhs_in  = gb.quicksum(y[r, j, i] for j in N if (j, i) in A)
                    z_out = gb.quicksum(z[(d, i, j)] for j in C if (i, j) in A and i != j and (self.base(i), self.base(j)) in Departures for d in Departures[(self.base(i), self.base(j))].keys())
                    z_in  = gb.quicksum(z[(d, j, i)] for j in C if (i, j) in A and i != j and (self.base(j), self.base(i)) in Departures for d in Departures[(self.base(j), self.base(i))].keys())
                    self.m.addConstr(lhs_out + z_out - lhs_in - z_in == fi_r[(r, i)], name=f"pass_bal_transfer[{r},{i}]")
        else: 
            if Cr:
                for r in R:
                    for i in Cr[r]:
                        self.m.addConstr(
                            gb.quicksum(y[r,i,j] for j in N if (i,j) in A) + gb.quicksum(z[(i,j)] for j in Cr[r] if (i,j) in A)
                            - gb.quicksum(y[r,j,i] for j in N if (j,i) in A) - gb.quicksum(z[(j,i)] for j in Cr[r]if (j,i) in A)
                            == fi_r[(r,i)],
                            name=f"pass_bal_transfer[{r},{i}]"
                        )

        # non-transfer nodes
        for r in R:
            for i in N:
                if i not in C:
                    self.m.addConstr(
                        gb.quicksum(y[r,i,j] for j in N if (i,j) in A) - gb.quicksum(y[r,j,i] for j in N if (j,i) in A) == fi_r[(r,i)],
                        name=f"pass_bal_nontransfer[{r},{i}]"
                    )

    def add_capacity_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, y = self.vars_['x'], self.vars_["v"], self.vars_["y"]
        qr, Q = self.params["qr"], self.params["Q"]
        if self.variable_substitution:
            for (i,j) in A:
                self.m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * v[i,j])

        else:
            for (i,j) in A:
                self.m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * gb.quicksum(x[k,i,j] for k in K))

    def add_time_consistency_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, z, T_node, T_veh, D_node = self.vars_["x"], self.vars_["v"], self.vars_["z"], self.vars_["T_node"], self.vars_["T_veh"], self.vars_["D_node"]
        tij, di, ei, li, ek, lk, Lbar, M, T, Cr = self.params["tij"], self.params["di"], self.params["ei"], self.params["li"], self.params["ek"], self.params["lk"], self.params['Lbar'], self.params["M"], self.params["T"], self.params['Cr']
        pair_pi_di = self.params["pair_pi_di"]
        zeroDepot_node, endDepot_node = zeroDepot, endDepot
        
        if self.variable_substitution:
            # depot to first node
            for k in K:
                for j in N:
                    if (zeroDepot_node,j) in A:
                        self.m.addConstr(T_node[j] >= T_veh[k] + tij[(self.base(zeroDepot_node), self.base(j))] - M*(1-x[k,zeroDepot_node,j]))

            # service times
            for (i,j) in A:
                self.m.addConstr(T_node[j] >= T_node[i] + di[self.base(i)] + tij[(self.base(i),self.base(j))] - M*(1-v[i,j]))

            # ride times
            for p, d in pair_pi_di.items():
                self.m.addConstr(T_node[d] - (T_node[p] + di[self.base(p)]) >= tij[(self.base(p), self.base(d))])
                self.m.addConstr(T_node[d] - (T_node[p] + di[self.base(p)]) <= Lbar[self.base(p)])

            # time windows
            for i in P+D:
                self.m.addConstr(T_node[i] >= ei[self.base(i)])
                self.m.addConstr(T_node[i] <= li[self.base(i)])

            # cumulative time
            for (i,j) in A:
                self.m.addConstr(D_node[j] >= D_node[i] + di[self.base(i)] + tij[(self.base(i),self.base(j))] - M*(1-v[i,j]))

            # max route duration
            for i in N:
                if (i,endDepot_node) in A:
                    self.m.addConstr(D_node[i] + di[self.base(i)] + tij[(self.base(i),self.base(endDepot_node))] <= T)

            # depot departure window
            for k in K:
                self.m.addConstr(T_veh[k] >= ek[k])
                self.m.addConstr(T_veh[k] <= lk[k])

            for r in R:
                if Cr:
                    for i in Cr[r]:
                        for j in Cr[r]:
                            if (i,j) in tij:
                                self.m.addConstr(
                                    T_node[j] >= T_node[i] + tij[self.base(i),self.base(j)] - M*(1 - z[(i,j)]),
                                    name=f"time_link_fixed[{r},{i},{j}]"
                                )
    
        else:
            # (9) Timing from depot to first served node on vehicle k
            for k in K:
                for i in N:
                    for j in N:
                        if (i,j) in tij and (zeroDepot_node,j) in tij:
                            self.m.addConstr(
                                T_node[j] >= T_veh[k] + tij[self.base(zeroDepot_node),self.base(j)] - M*(1 - x[k,zeroDepot_node,j]),
                                name=f"time_dep_first[{k},{j}]"
                            )

            # (10) Timing between consecutive service nodes (DRT)
            for i in N:
                for j in N:
                    if (self.base(i),self.base(j)) in tij:
                        self.m.addConstr(
                            T_node[j] >= T_node[i] + di[self.base(i)] + tij[self.base(i),self.base(j)] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                            name=f"time_link_drt[{i},{j}]"
                        )

            # (11) Timing along fixed-route between transfer nodes for request r
            for r in R:
                if Cr:
                    for i in Cr[r]:
                        for j in Cr[r]:
                            if (i,j) in A:
                                self.m.addConstr(
                                    T_node[j] >= T_node[i] + tij[self.base(i),self.base(j)] - M*(1 - z[(i,j)]),
                                    name=f"time_link_fixed[{r},{i},{j}]"
                                )

            # (12) Time windows at pickup and dropoff nodes
            ### for i in (P + D): ### No time windows for dropoffs in this example
            for i in P + D:
                self.m.addConstr(T_node[i] >= ei[self.base(i)], name=f"tw_lo[{i}]")
                self.m.addConstr(T_node[i] <= li[self.base(i)], name=f"tw_hi[{i}]")

            # (13) Ride-time constraints: t_{i,n+i} <= B_{n+i} - (B_i + d_i) <= Lbar_i  for pickups i
            for i in P:
                drop = pair_pi_di[i]  # dropoff node for pickup i (usually n+i)
                if (i,drop) in A:
                    self.m.addConstr(T_node[drop] - (T_node[i] + di[self.base(i)]) >= tij[self.base(i),self.base(drop)], name=f"ride_lo[{i}]")
                    self.m.addConstr(T_node[drop] - (T_node[i] + di[self.base(i)]) <= Lbar[self.base(i)],   name=f"ride_hi[{i}]")

            # (14) Cumulative time variable D (mirrors (10) without B, used for duration in (15))
            for i in N:
                for j in N:
                    if (i,j) in A:
                        self.m.addConstr(
                            D_node[j] >= D_node[i] + di[self.base(i)] + tij[self.base(i),self.base(j)] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                            name=f"D_link[{i},{j}]"
                        )

            # (15) Vehicle maximum route duration T (ensure from node i to end depot fits in T)
            for i in N:
                if (i,endDepot_node) in A:
                    self.m.addConstr(D_node[i] + di[self.base(i)] + tij[self.base(i),self.base(endDepot_node)] <= T, name=f"duration[{i}]")

            # (16) Depot departure window for each vehicle
            for k in K:
                self.m.addConstr(T_veh[k] >= ek[k], name=f"veh_dep_lo[{k}]")
                self.m.addConstr(T_veh[k] <= lk[k], name=f"veh_dep_hi[{k}]")

    def add_variable_substitution_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        v, x = self.vars_["v"], self.vars_["x"]
        zeroDepot_node = zeroDepot

        for (i,j) in A:
            if i != zeroDepot_node:
                self.m.addConstr(v[i,j] == gb.quicksum(x[k,i,j] for k in K))

    def add_transfer_node_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, z = self.vars_['x'], self.vars_["v"], self.vars_["z"]
        Cr, Departures = self.params['Cr'], self.params["Departures"]
        
        if self.variable_substitution:
            if self.timetabled_departures and Departures != None:
                for i in C:
                    self.m.addConstr(
                        gb.quicksum(v[j,i] for j in N if (j,i) in A) ==
                        gb.quicksum(z[(d,i,j)] + z[(d,j,i)] for j in C if (i,j) in A and i!=j and (self.base(i),self.base(j)) in Departures for d in Departures[(self.base(i),self.base(j))].keys())
                    )
            else:
                if Cr:
                    for r in R:
                        for i in Cr[r]:
                            self.m.addConstr(
                                gb.quicksum(v[j,i] for j in N if (j,i) in A) ==
                                gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) in A)
                            )
                    for r in R:
                        self.m.addConstr(gb.quicksum(v[i,j] for i in Cr[r] for j in N if (i,j) in A) <= 2)
        else:
            if Cr:
                for r in R:
                    for i in Cr[r]:
                        self.m.addConstr(
                            gb.quicksum(gb.quicksum(x[k,j,i] for k in K) for j in N if (j,i) in A) ==
                            gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) in A)
                        )
            if Cr:
                for r in R:
                    self.m.addConstr(gb.quicksum(gb.quicksum(x[k,i,j] for k in K) for i in Cr[r] for j in N if (i,j) in A) <= 2)    
        
    def add_battery_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        x, v, B_node, E_node = self.vars_['x'], self.vars_["v"], self.vars_["B_node"], self.vars_["E_node"]
        
        tij = self.params["tij"]
        beta, alpha, C_bat, gamma = self.params["beta"], self.params["alpha"], self.params["C_bat"], self.params["gamma"]

        if self.variable_substitution:
            for (i, j) in A:
                if i not in F and j not in F:
                    self.m.addConstr(B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] - C_bat * (1 - v[i, j]) <= 0)
                    self.m.addConstr(B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] + C_bat * (1 - v[i, j]) >= 0)
        else:
            for (i, j) in A:
                if i not in F and j not in F:
                    self.m.addConstr(B_node[j] - B_node[i] + beta * tij[(self.base(i), self.base(j))] - C_bat * (1 - gb.quicksum(x[k,i,j] for k in K)) <= 0)

        # battery capacity
        for s in F:
            self.m.addConstr(C_bat - (B_node[s] + alpha * E_node[s]) >= 0)

        # min SOC
        for i in N:
            self.m.addConstr(B_node[i] >= gamma * C_bat)

        # init SOC
        self.m.addConstr(B_node[zeroDepot] == C_bat)

    def add_scheduled_PT_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        z, T_node = self.vars_["z"], self.vars_["T_node"]
        M, Departures, tij = self.params["M"], self.params["Departures"], self.params['tij']

        for i in N:
            if i not in C:
                self.m.addConstr(
                    M*(1 - gb.quicksum(z[(d,i,j)] for j in C if (i, j) in A and i != j and (self.base(j), self.base(i)) in Departures for d in Departures[(self.base(i),self.base(j))].keys())) +
                    gb.quicksum(Departures[(self.base(i),self.base(j))][d] * z[(d,i,j)] for j in C if (i,j) in A and i!=j and (self.base(i),self.base(j)) in Departures for d in Departures[(self.base(i),self.base(j))].keys())
                    - T_node[i] >= 0
                )
        for i in C:
            expr = gb.quicksum(
                (Departures[(self.base(i), self.base(j))][d] + tij[self.base(i), self.base(j)]) * z[d, i, j]
                for j in C
                if (i, j) in A
                and (self.base(i), self.base(j)) in Departures
                for d in Departures[(self.base(i), self.base(j))].keys()
            )
            self.m.addConstr(T_node[i] - expr >= 0, name=f"PT_departure_after_arrival_{i}")

    def add_scheduled_PT_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        z, a, T_node, y = self.vars_["z"], self.vars_["a"], self.vars_["T_node"], self.vars_["y"]
        M, tij, Departures = self.params["M"], self.params['tij'], self.params["Departures"]

        # Constraint : PT departure after service start at transfer node
        for i in C:
            expr = (M + gb.quicksum((Departures[self.base(i), self.base(j)][d] - M) * z[d, i, j] for j in C for d in Departures[self.base(i), self.base(j)]) - T_node[i])
            self.m.addConstr(expr >= 0, name=f"PT_departure_after_service_{i}")

        # Constraint : Service at arrival transfer node after PT arrival
        for j in C:
            expr = (T_node[j] - gb.quicksum((Departures[self.base(i), self.base(j)][d] + tij[self.base(i), self.base(j)]) * z[d, i, j] for i in C for d in Departures[self.base(i), self.base(j)]))
            self.m.addConstr(expr >= 0, name=f"PT_arrival_before_service_{j}")

        # Constraint : If node used in PT trip → mark as visited
        for i in C:
            expr_lhs = M * a[i]
            expr_rhs = (
                gb.quicksum(z[d, i, j] for j in C for d in Departures[self.base(i), self.base(j)]) +
                gb.quicksum(z[d, j, i] for j in C for d in Departures[self.base(j), self.base(i)])
            )
            self.m.addConstr(expr_lhs >= expr_rhs, name=f"PT_node_visit_upper_{i}")

        # Constraint : Node marked as visited only if used in PT trip
        for i in C:
            expr_rhs = (
                gb.quicksum(z[d, i, j] for j in C for d in Departures[self.base(i), self.base(j)]) +
                gb.quicksum(z[d, j, i] for j in C for d in Departures[self.base(j), self.base(i)])
            )
            self.m.addConstr(a[i] <= expr_rhs, name=f"PT_node_visit_lower_{i}")

        # Constraint : Linearised constraint such that Ti + tij <= Tj when going from i to j for request r
        for j1 in C:
            for j2 in C:
                if (self.base(j1), self.base(j2)) in Departures:
                    for r in R:
                        lhs = T_node[j2]
                        rhs = (
                            T_node[j1] + tij[self.base(j1), self.base(j2)]
                            - M * (
                                3
                                - gb.quicksum(y[r, i, j1] for i in N if i not in C)
                                - gb.quicksum(y[r, j2, i] for i in N if i not in C)
                                - gb.quicksum(z[d, j1, j2] for d in Departures[self.base(j1), self.base(j2)].keys())
                            )
                        )
                        self.m.addConstr(lhs >= rhs, name=f"PT_DARP_sync_{r}_{j1}_{j2}")

                        
        # for r in R:
        #     for i1 in N:
        #         for i2 in N:
        #             for j1 in C:
        #                 for j2 in C:
        #                     if y[r,i1,j1].X * y[r,i2,j2].X > 1e-3:
        #                         self.m.addConstr(T_node[j2] - T_node[j1] - tij[self.base(j1), self.base(j2)] >= 0)

        # for i in C:
        #     for j in C:
        #         for d in Departures[self.base(i) , self.base(j)]:
        #             if z[d,i,j].x > 0.5 : 
        #                 self.m.addConstr(T_node[j2] - T_node[j1] - tij[self.base(j1), self.base(j2)] >= 0)




    def add_artificial_node_constraints(self):
        nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
        M = self.params['M']
        z, a, T_node = self.vars_["z"], self.vars_["a"], self.vars_['T_node']
        for (i,m_) in N:
            if m_ < max(n for (j,n) in N if j == i):
                self.m.addConstr(a[i,m_] - a[i,m_+1] >= 0)
                self.m.addConstr(T_node[i,m_ + 1] - T_node[i,m_] >= 0) # + M * (2 - a[i,m_] - a[i,m_+1]) <= 0, name = "Time Logic for artificial nodes")




    # def add_MoPS_constraints(self):
    #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = self.sets["nodes"], self.sets["N"], self.sets["P"], self.sets["D"], self.sets["C"], self.sets["F"], self.sets["R"], self.sets["K"], self.sets["zeroDepot"], self.sets["endDepot"], self.sets["A"]
    #     v = self.vars_['v']
    #     pair_pi_di = self.params['pair_pi_di']
    #     for i in P_M + D_M:
    #         self.m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in A) <= 1, name="Each MoPS entered less than 1.5 (exit)")
    #         self.m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) <= 1, name="Each MoPS entered less than 1.5 (entry)")

    #     for i in P_M:
    #         self.m.addConstr(gb.quicksum(v[i,j] for j in N in (i,j) in A) - gb.quicksum(v[pair_pi_di(i), j] for j in N if (pair_pi_di(i), j) in A) == 0, name="Every Picked-up MoPS request must be dropped off")

