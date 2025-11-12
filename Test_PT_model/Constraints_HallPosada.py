##### Imports #####
import gurobipy as gb

class DARPConstraintBuilder_HallPosada:
    """
    Reimplementation of Model 2 from:
        Häll & Posada (2017)
    Adapted for:
        - Homogeneous fleet (single class)
        - Single resource
        - Timetabled PT arcs
    """

    def __init__(self, m, vars_, sets, params):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params

    def base(self, x):
        """Return base node id (flat int) from (i,m) or i."""
        return x[0] if isinstance(x, tuple) else x

    # --------------------------------------------------------------------------
    #  VARIABLES
    # --------------------------------------------------------------------------
    def define_variables(self):
        """Define all decision variables for Model 2."""
        N, NG, NP, ND, A, K, R = (
            self.sets["N"], self.sets["NG"], self.sets["NP"], self.sets["ND"],
            self.sets["A"], self.sets["K"], self.sets["R"]
        )
        Dij = self.params["Dij"]

        # Vehicle routing arcs
        x = self.m.addVars(A, K, vtype=gb.GRB.BINARY, name="x[imjnk]")
        # Passenger flows
        # Passenger flow arcs (exclude transfer→transfer arcs)
        NG_nodes = {i for (i, m) in NG}
        DAR_arcs = [
            (i, j) for (i, j) in A
            if not (self.base(i) in NG_nodes and self.base(j) in NG_nodes)
        ]
        y = self.m.addVars(R, A, vtype=gb.GRB.BINARY, name="y[imjnr]")

        # PT arcs (departures)
        z = {}
        for (i, j) in Dij:
            for d in Dij[(i, j)]:
                for r in R:
                    z[(i, j, d, r)] = self.m.addVar(vtype=gb.GRB.BINARY,
                                                    name=f"z[{i},{j},{d},{r}]")
        # Walking arcs (unused but kept)
        w = self.m.addVars(R, A, vtype=gb.GRB.BINARY, name="w[imjnr]")
        # Transfer activation
        a = self.m.addVars(NG, vtype=gb.GRB.BINARY, name="a[im]")
        # Service start time
        t = self.m.addVars(N, vtype=gb.GRB.CONTINUOUS, lb=0, name="t[im]")

        self.vars_ = {"x": x, "y": y, "z": z, "w": w, "a": a, "t": t}
        return self.vars_

    # --------------------------------------------------------------------------
    #  OBJECTIVE (Eq. 19)
    # --------------------------------------------------------------------------
    def set_objective(self):
        x, z = self.vars_["x"], self.vars_["z"]
        A, K = self.sets["A"], self.sets["K"]
        Cij, Cij_dr = self.params["Cij"], self.params["Cij_dr"]

        # Vehicle cost
        obj = gb.quicksum(
            Cij[self.base(i), self.base(j)] * x[i, j, k]
            for (i, j) in A for k in K
            if (self.base(i), self.base(j)) in Cij
        )
        # PT cost
        obj += gb.quicksum(
            Cij_dr[self.base(i), self.base(j), d, r] * z[i, j, d, r]
            for (i, j, d, r) in z.keys()
            if (self.base(i), self.base(j), d, r) in Cij_dr
        )

        self.m.setObjective(obj, gb.GRB.MINIMIZE)

    # --------------------------------------------------------------------------
    #  VEHICLE CONSTRAINTS (20–23)
    # --------------------------------------------------------------------------
    def add_vehicle_constraints(self):
        x = self.vars_["x"]
        A, NP, ND, N, K = self.sets["A"], self.sets["NP"], self.sets["ND"], self.sets["N"], self.sets["K"]
        startDepot, endDepot = self.sets["zeroDepot"], self.sets["endDepot"]

        # (20) each pickup/dropoff visited once
        for i in NP + ND:
            self.m.addConstr(
                gb.quicksum(x[ii, jj, k] for (ii, jj) in A if ii == i for k in K) == 1,
                name=f"(20)_visit_once_{i}"
            )

        # (21) each vehicle leaves depot once
        for k in K:
            self.m.addConstr(
                gb.quicksum(x[startDepot, j, k] for (i, j) in A if i == startDepot) == 1,
                name=f"(21)_depot_start_{k}"
            )

        # (22) flow conservation
        for (i, m) in N:
            if (i, m) not in [startDepot, endDepot]:
                for k in K:
                    inflow = gb.quicksum(x[j, (i, m), k] for (j, jj) in A if jj == (i, m))
                    outflow = gb.quicksum(x[(i, m), j, k] for (ii, j) in A if ii == (i, m))
                    self.m.addConstr(inflow - outflow == 0, name=f"(22)_flow_cons_{i}_{k}")

        # (23) each vehicle ends at depot
        for k in K:
            self.m.addConstr(
                gb.quicksum(x[i, endDepot, k] for (i, j) in A if j == endDepot) == 1,
                name=f"(23)_depot_end_{k}"
            )

    # --------------------------------------------------------------------------
    #  PASSENGER BALANCE (24–26)
    # --------------------------------------------------------------------------
    def add_passenger_balance_constraints(self):
        """
        Passenger flow balance constraints (Eqs. 24–25).
        Ensures each passenger request r has balanced inflow/outflow
        at pickup, drop-off, and transfer nodes.
        """
        y, z = self.vars_["y"], self.vars_["z"]
        A, NG, NP, ND, R = (
            self.sets["A"], self.sets["NG"],
            self.sets["NP"], self.sets["ND"], self.sets["R"]
        )
        Fir = self.params["Fir"]

        # (24) Pickup / drop-off balance: outflow − inflow = F_ir
        for r in R:
            for i in NP + ND:
                inflow = gb.quicksum(y[r, j, i] for (j, jj) in A if jj == i)
                outflow = gb.quicksum(y[r, i, j] for (ii, j) in A if ii == i)
                rhs = Fir[r, i]
                self.m.addConstr(outflow - inflow == rhs,
                                 name=f"(24)_balance_{r}_{i}")

        # (25) Transfer-node balance: total inflow (y + z) = total outflow (y + z)
        for r in R:
            for (i, m) in NG:
                # --- DAR arcs ---
                inflow_y = gb.quicksum(y[r, j, jj] for (j, jj) in A if jj == (i, m))
                outflow_y = gb.quicksum(y[r, ii, j] for (ii, j) in A if ii == (i, m))

                # --- PT arcs (flat indices) ---
                inflow_z = gb.quicksum(
                    z[jj, self.base(i), d, r]
                    for (jj, ii, d, rr) in z.keys()
                    if ii == self.base(i) and rr == r
                )
                outflow_z = gb.quicksum(
                    z[self.base(i), jj, d, r]
                    for (ii, jj, d, rr) in z.keys()
                    if ii == self.base(i) and rr == r
                )

                self.m.addConstr(
                    (outflow_y + outflow_z) - (inflow_y + inflow_z) == Fir[r, (i,m)],
                    name=f"(25)_transfer_balance_{r}_{i}_{m}"
                )


    # --------------------------------------------------------------------------
    #  CAPACITY (27)
    # --------------------------------------------------------------------------
    def add_capacity_constraints(self):
        x, y = self.vars_["x"], self.vars_["y"]
        A, K, R = self.sets["A"], self.sets["K"], self.sets["R"]
        Lr, Q = self.params["Lr"], self.params["Q"]

        for (i, j) in [(i, j) for (i, j, k) in x.keys()]:
            self.m.addConstr(
                gb.quicksum(Lr[self.base(r)] * y[r, i, j] for r in R) <=
                Q * gb.quicksum(x[i, j, k] for k in K),
                name=f"(27)_capacity_{i}_{j}"
            )

    # --------------------------------------------------------------------------
    #  TIME CONSTRAINTS (28–36)
    # --------------------------------------------------------------------------
    def add_time_constraints(self):
        x, z, t = self.vars_["x"], self.vars_["z"], self.vars_["t"]
        A, K, NG, NP, ND = self.sets["A"], self.sets["K"], self.sets["NG"], self.sets["NP"], self.sets["ND"]
        M = self.params["M"]
        rbar = self.params.get("rbar", 4)
        T0j, Tij, TD, TA, Bi, Ti_e, Ti_l = (
            self.params["T0j"], self.params["Tij"], self.params["TD"], self.params["TA"],
            self.params["Bi"], self.params["Ti_e"], self.params["Ti_l"]
        )
        startDepot = self.sets["zeroDepot"]

        # (28) depot to first node
        for (j, n) in [node for node in t.keys() if node != startDepot]:
            if self.base(j) in T0j:
                self.m.addConstr(
                    t[j, n] - T0j[self.base(j)] +
                    M * (1 - gb.quicksum(x[startDepot, j, k] for k in K)) >= 0,
                    name=f"(28)_depot_to_first_{j}"
                )

        # (29) service time consistency
        for (i, j, k) in x.keys():
            if (self.base(i), self.base(j)) in Tij:
                self.m.addConstr(
                    t[j] - t[i] - Tij[self.base(i), self.base(j)] +
                    M * (1 - x[i, j, k]) >= 0,
                    name=f"(29)_time_progress_{i}_{j}_{k}"
                )

        # (32) PT departure after service start
        for (i, m) in NG:
            for r in self.sets["R"]:
                lhs = (M +
                       gb.quicksum(
                           (TD[self.base(i), self.base(j), d] - M) * z[i, j, d, r]
                           for (i, j, d, rr) in z.keys() if i == self.base(i) and rr == r
                           if (self.base(i), self.base(j), d) in TD
                       )
                       - t[i, m])
                self.m.addConstr(lhs >= 0, name=f"(32)_PT_depart_after_{i}_{r}")

        # (33) PT arrival before service start
        for (j, n) in NG:
            for r in self.sets["R"]:
                lhs = (t[j, n] -
                       gb.quicksum(
                           (TA[self.base(i), self.base(j), d]) * z[i, j, d, r]
                           for (i, j, d, rr) in z.keys() if j == self.base(j) and rr == r
                           if (self.base(i), self.base(j), d) in TA
                       ))
                self.m.addConstr(lhs >= 0, name=f"(33)_PT_arrival_before_{j}_{r}")

        # (34) transfer time order
        for (i, m) in NG:
            if (i, m + 1) in NG:
                self.m.addConstr(t[i, m] - t[i, m + 1] <= 0,
                                 name=f"(34)_transfer_time_order_{i}_{m}")

        # (35) time windows
        for (i, m) in NP + ND:
            base_i = self.base(i)
            if base_i in Ti_e and base_i in Ti_l:
                self.m.addConstr(t[i, m] >= Ti_e[base_i], name=f"(35a)_timewin_lower_{i}")
                self.m.addConstr(t[i, m] <= Ti_l[base_i], name=f"(35b)_timewin_upper_{i}")

        # (36) ride time
        for (i, m) in NP:
            drop_id = (i + rbar, 1)
            if drop_id in t.keys() and self.base(i) in Bi:
                self.m.addConstr(
                    t[drop_id] - t[i, m] <= Bi[self.base(i)],
                    name=f"(36)_ride_time_{i}"
                )

    # --------------------------------------------------------------------------
    #  TRANSFER NODE CONSTRAINTS (37–39)
    # --------------------------------------------------------------------------    
    def add_transfer_node_constraints(self):
        """
        (37)–(39) Transfer node activation consistency.
        Ensures that transfer node activation variable a[i,m]
        is consistent with PT flow usage z[i,j,d,r].
        """
        z, a = self.vars_["z"], self.vars_["a"]
        NG = self.sets["NG"]
        M = self.params["M"]

        for (i, m) in NG:
            # total PT flow (sum of all z arcs entering or leaving base(i))
            total_z = gb.quicksum(
                z[ii, jj, d, r]
                for (ii, jj, d, r) in z.keys()
                if ii == self.base(i) or jj == self.base(i)
            )

            # (37) Upper bound: if transfer active, it can support PT flow
            self.m.addConstr(
                M * a[i, m] >= total_z,
                name=f"(37)_node_visit_upper_{i}_{m}"
            )

            # (38) Lower bound: if any PT flow exists, a[i,m] must be active
            # Linear safe form: a[i,m] ≤ total_z (scaled if necessary)
            self.m.addConstr(
                a[i, m] <= total_z + 1e-6,
                name=f"(38)_node_visit_lower_{i}_{m}"
            )

            # (39) Monotonic activation order (higher m ⇒ not less active)
            if (i, m + 1) in NG:
                self.m.addConstr(
                    a[i, m + 1] - a[i, m] <= 0,
                    name=f"(39)_activation_order_{i}_{m}"
                )


    def add_strengthening_constraints(self):
        """
        Adds strengthening constraints (47)–(52) from Häll & Posada (2017)
        to tighten the LP relaxation.
        Subtour elimination is NOT included.
        """
        x, y, z, a, t = (
            self.vars_["x"], self.vars_["y"],
            self.vars_["z"], self.vars_["a"], self.vars_["t"]
        )
        sets, params = self.sets, self.params
        A, NP, ND, NG, R, K = sets["A"], sets["NP"], sets["ND"], sets["NG"], sets["R"], sets["K"]
        M = params["M"]

        # ------------------------------------------------------------------
        # (47) Link between passenger flow and vehicle arcs
        # ------------------------------------------------------------------
        #  y_{r,ij} <= ∑_k x_{ijk}   ∀ (i,j)∈A, ∀r∈R
        for r in R:
            for (i, j) in A:
                self.m.addConstr(
                    y[r, i, j] <= gb.quicksum(x[i, j, k] for k in K),
                    name=f"(47)_y_le_x_{r}_{i}_{j}"
                )

        # ------------------------------------------------------------------
        # (48) Pickup before drop-off timing reinforcement
        # ------------------------------------------------------------------
        #  t_d >= t_p + s_p  (service sequence for each request)
        rbar = params.get("rbar", 4)
        for (i, m) in NP:
            drop = (i + rbar, 1)
            if drop in t.keys():
                self.m.addConstr(
                    t[drop] >= t[i, m],
                    name=f"(48)_pickup_before_drop_{i}"
                )

        # ------------------------------------------------------------------
        # (49) Time-window propagation through active arcs
        # ------------------------------------------------------------------
        #  t_j >= t_i + (Tij_ij - M*(1 - x_ijk))
        Tij = params["Tij"]
        for (i, j, k) in x.keys():
            if (self.base(i), self.base(j)) in Tij:
                self.m.addConstr(
                    t[j] >= t[i] + Tij[self.base(i), self.base(j)] - M * (1 - x[i, j, k]),
                    name=f"(49)_time_propagation_{i}_{j}_{k}"
                )

        # ------------------------------------------------------------------
        # (50) Transfer activation timing strengthening
        # ------------------------------------------------------------------
        #  t_{i,m+1} >= t_{i,m}  when a[i,m+1] active
        for (i, m) in NG:
            if (i, m + 1) in NG:
                self.m.addConstr(
                    t[i, m + 1] >= t[i, m] - M * (1 - a[i, m + 1]),
                    name=f"(50)_transfer_time_strength_{i}_{m}"
                )

        # ------------------------------------------------------------------
        # (51) No-service without activation (tightening for a)
        # ------------------------------------------------------------------
        #  t_{i,m} <= M * a_{i,m}
        for (i, m) in NG:
            self.m.addConstr(
                t[i, m] <= M * a[i, m],
                name=f"(51)_no_service_without_activation_{i}_{m}"
            )

        # ------------------------------------------------------------------
        # (52) Vehicle arrival consistency tightening
        # ------------------------------------------------------------------
        #  If x_{ijk} active → enforce earliest arrival window
        Ti_e = params.get("Ti_e", {})
        for (i, j, k) in x.keys():
            base_i = self.base(i)
            if base_i in Ti_e:
                self.m.addConstr(
                    t[i] >= Ti_e[base_i] - M * (1 - x[i, j, k]),
                    name=f"(52)_arrival_window_{i}_{j}_{k}"
                )

        print("✅ Strengthening constraints (47–52) added.")

    def add_extra_constraints(self):
        x, y = self.vars_["x"], self.vars_["y"]
        A, R, K = self.sets["A"], self.sets["R"], self.sets["K"]
        for r in R:
            for (i, j) in A:
                self.m.addConstr(
                    y[r, i, j] <= gb.quicksum(x[i, j, k] for k in K),
                    name=f"(PV)_y_le_x_r{r}_i{self.base(i)}_j{self.base(j)}"
                )



    # --------------------------------------------------------------------------
    #  FINALIZE
    # --------------------------------------------------------------------------
    def finalize(self):
        self.m.update()
        self.m.write("HallPosada_Model2.lp")
        print("✅ Model 2 (Häll & Posada 2017) constraints successfully built.")
