##### Imports #####

import gurobipy as gb

import Graph_Builder_3 as g

def build_idarp_model_Enhanced_5_EV(Time_Limit):

    data = g.data   # grab from Graph_builder

    # Unpack sets
    N, P, D, C, R, K, Cr, F = data['N'], data['P'], data['D'], data['C'], data['R'], data['K'], data['Cr'], data['F']
    Q, T, cij, tij, qr, di, ei, li, ek, lk = data['Q'], data['T'], data['cij'], data['tij'], data['qr'], data['di'], data['ei'], data['li'], data['ek'], data['lk']
    fi_r, Lbar, endDepot, zeroDepot, pair_pi_di, M, n, n_K = data['fi_r'], data['Lbar'], data['endDepot'], data['zeroDepot'], data['pair_pi_di'], data['M'], data['n'], data['n_K']
    C_bat, alpha, beta, gamma = data['C_bat'], data['alpha'], data['beta'], data['gamma']

    # Build arc list from tij keys
    A = [(i, j) for (i, j) in tij.keys()]
    Aset = set(A)

    m = gb.Model("IDARP_Enhanced_5_EV")
    m.setParam('TimeLimit', Time_Limit)

    # Decision variables
    # x_{ij}^k ∈ {0,1}
    x = m.addVars(K, A, vtype=gb.GRB.BINARY, name="x")

    # v_{ij} ∈ {0,1}
    v = m.addVars(A, vtype=gb.GRB.BINARY, name="v")

    # y_{ij}^r ≥ 0  (implicitly integral in feasible solutions)
    y = m.addVars(R, A, vtype=gb.GRB.CONTINUOUS, lb=0, name="y")

    # z_{ij} for (i,j) in C_r × C_r — we index it as z[r,i,j] to avoid ambiguity
    z = {}
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                z[(i,j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")

    # B_i: start of service at node i
    T_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")

    # B_k: time vehicle k leaves the depot
    T_veh  = m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")

    # D_i: time elapsed from vehicle departure until arrival at node i
    D_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

    # B_i: Battery state at node i.
    B_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")

    # E_i: Time spent recharging at station i in F
    E_node = m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")

    # Objective (1): min sum c_ij * sum_k x_ij^k
    m.setObjective(gb.quicksum(cij[0,j] * x[k,0,j] for k in K for j in N if (0,j) in Aset) + gb.quicksum(cij[i,j] * v[i,j] for (i,j) in Aset if i!=0), gb.GRB.MINIMIZE)

    # ---------- Constraints ----------

    # (2) Every pickup/dropoff node visited exactly once by some vehicle
    for i in P + D:
        m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in Aset) == 1, name=f"visit_once[{i}]")
        m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in Aset) == 1, name=f"visit_once_rev[{i}]")

    # (2 - annex) Every transfer node visited at most once by some vehicle
    for i in C:
        m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in Aset) <= 1, name=f"visit_once_transfer[{i}]")
        m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in Aset) <= 1, name=f"visit_once_transfer_rev[{i}]")

    # (3) Each vehicle leaves the origin depot exactly once (link x and v)
    for k in K:
        m.addConstr(gb.quicksum(x[k,zeroDepot,j] for j in N if (zeroDepot,j) in Aset) <= 1, name=f"start_from_depot[{k}]")
    for j in N:
        if (zeroDepot, j) in Aset:
            m.addConstr(gb.quicksum(x[k, zeroDepot, j] for k in K) == v[zeroDepot, j],
                        name=f"link_depot_x_v_to_v[0,{j}]")


    # (4) Vehicle flow conservation at pickups, dropoffs, and transfers
    for k in K:
        for i in (P + D + C + F):
            m.addConstr(gb.quicksum(x[k,j,i] for j in N if (j,i) in Aset) - gb.quicksum(x[k,i,j] for j in N if (i,j) in Aset) == 0,
                        name=f"flow_cons_x[{k},{i}]")
    for i in (P + D + C + F):
            m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in Aset) - gb.quicksum(v[i,j] for j in N if (i,j) in Aset) == 0,
                        name=f"flow_cons_v[{i}]")

    # (5) Each vehicle ends at the destination depot exactly once
    m.addConstr(gb.quicksum(v[i,endDepot] for i in N if (i,endDepot) in Aset) == n_K, name=f"end_at_depot[{k}]")
    for k in K:
        m.addConstr(gb.quicksum(x[k,i,endDepot] for i in N if (i,endDepot) in Aset) == 1, name=f"veh_end_at_depot[{k}]")

    # (6) Passenger balance at transfer nodes (can be via DRT or fixed route)
    for r in R:
        for i in Cr[r]:
            m.addConstr(
                gb.quicksum(y[r,i,j] for j in N if (i,j) in Aset) + gb.quicksum(z[(i,j)] for j in Cr[r] if (i,j) in Aset)
                - gb.quicksum(y[r,j,i] for j in N if (j,i) in Aset) - gb.quicksum(z[(j,i)] for j in Cr[r]if (j,i) in Aset)
                == fi_r[(r,i)],
                name=f"pass_bal_transfer[{r},{i}]"
            )

    # (7) Passenger balance at non-transfer nodes (only DAR arcs)
    for r in R:
        for i in N:
            if i not in Cr[r]:
                m.addConstr(
                    gb.quicksum(y[r,i,j] for j in N if (i,j) in Aset) - gb.quicksum(y[r,j,i] for j in N if (j,i) in Aset) == fi_r[(r,i)],
                    name=f"pass_bal_nontransfer[{r},{i}]"
                )

    # (8) Capacity on vehicle arcs
    for (i,j) in A:
        m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * v[i,j],
                    name=f"capacity[{i},{j}]")
        m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * gb.quicksum(x[k,i,j] for k in K),
                    name=f"capacity_x[{i},{j}]")

    # (9) Timing from depot to first served node on vehicle k
    for k in K:
        for j in N:
            if (zeroDepot,j) in Aset:
                m.addConstr(
                    T_node[j] >= T_veh[k] + tij[zeroDepot,j] - M*(1 - x[k,zeroDepot,j]),
                    name=f"time_dep_first[{k},{j}]"
                )

    # (10) Timing between consecutive service nodes (DRT)
    for (i,j) in A:
        m.addConstr(
            T_node[j] >= T_node[i] + di[i] + tij[i,j] - M*(1 - v[i,j]),
            name=f"time_link_drt[{i},{j}]")
        m.addConstr(
            T_node[j] >= T_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
            name=f"time_link_drt_x[{i},{j}]")

    # (11) Timing along fixed-route between transfer nodes for request r
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                if (i,j) in Aset:
                    m.addConstr(
                        T_node[j] >= T_node[i] + tij[i,j] - M*(1 - z[(i,j)]),
                        name=f"time_link_fixed[{r},{i},{j}]"
                    )

    # (12) Time windows at pickup and dropoff nodes
    ### for i in (P + D): ### No time windows for dropoffs in this example
    for i in P + D:
        m.addConstr(T_node[i] >= ei[i], name=f"tw_lo[{i}]")
        m.addConstr(T_node[i] <= li[i], name=f"tw_hi[{i}]")

    # (13) Ride-time constraints: t_{i,n+i} <= B_{n+i} - (B_i + d_i) <= Lbar_i  for pickups i
    for i in P:
        drop = pair_pi_di[i]  # dropoff node for pickup i (usually n+i)
        if (i,drop) in Aset:
            m.addConstr(T_node[drop] - (T_node[i] + di[i]) >= tij[i,drop], name=f"ride_lo[{i}]")
            m.addConstr(T_node[drop] - (T_node[i] + di[i]) <= Lbar[i],   name=f"ride_hi[{i}]")

    # (14) Cumulative time variable D (mirrors (10) without B, used for duration in (15))
    for (i,j) in A:
        m.addConstr(
            D_node[j] >= D_node[i] + di[i] + tij[i,j] - M*(1 - v[i,j]),
            name=f"D_link[{i},{j}]"
        )
        m.addConstr(
            D_node[j] >= D_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
            name=f"D_link_x[{i},{j}]")

    # (15) Vehicle maximum route duration T (ensure from node i to end depot fits in T)
    for i in N:
        if (i,endDepot) in Aset:
            m.addConstr(D_node[i] + di[i] + tij[i,endDepot] <= T, name=f"duration[{i}]")

    # (16) Depot departure window for each vehicle
    for k in K:
        m.addConstr(T_veh[k] >= ek[k], name=f"veh_dep_lo[{k}]")
        m.addConstr(T_veh[k] <= lk[k], name=f"veh_dep_hi[{k}]")

    ##### Variable Substitution Constraints #####
    # (17 - Added) One and only one vehicle can travel on arc (i,j)
    for (i,j) in A:
        if i != 0:
            m.addConstr(v[i,j] == gb.quicksum(x[k,i,j] for k in K), name=f"link_v_x[{i},{j}]")

    # (18 - added) Only one vehicle can go from depot to j
    for j in P+D+C:
        if (0,j) in Aset:
            m.addConstr(gb.quicksum(x[k,0,j] for k in K) == v[0,j], name=f"no_direct_depot[{j}]")

    ##### Subtour elimination constraints #####
    # (19 -added) Only one vehicle can travel between a pair of nodes and only in one direction
    for i in N:
        for j in N:
            if i != j and (i,j) in Aset and (j,i) in Aset:
                m.addConstr(v[i,j] + v[j,i] <= 1, name=f"subtour_elim[{i},{j}]")
    
    ##### Transfer node strengthening constraints #####
    # (20 - added) A maximum of two transfer nodes can be visited per request
    for r in R:
        m.addConstr(gb.quicksum(v[i,j] for i in Cr[r] for j in N if (i,j) in Aset) <= 2, name=f"max_transfers[{r}]")

    # (21 - added) If a vehicle visits a transfer node, it must use or have used the fixed-route arc for that request
    for r in R:
        for i in Cr[r]:
            m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in Aset) == 
                        gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) in z),
                        name=f"transfer_sync[{r},{i}]")
            
    # ##### EV constraints #####
    # # Battery depletion rate
    # for i in N:
    #     for j in N:
    #         if i not in F:
    #             if j not in F:
    #                 if (i,j) in Aset:
    #                     m.addConstr(B_node[j] - B_node[i] + beta*tij[i,j] - C_bat*(1 - v[i,j]) <= 0, name=f"battery_depletion_ub")
    #                     m.addConstr(B_node[j] - B_node[i] + beta*tij[i,j] + C_bat*(1 - v[i,j]) >= 0, name=f"battery_depletion_lb")

    # # Battery Charging rate 
    # for j in N:
    #     for s in F:
    #         if j not in F:
    #             if (s,j) in Aset:
    #                 if j != s:
    #                     m.addConstr(B_node[j] - B_node[s] + beta*tij[s,j] - alpha*E_node[s] - C_bat*(1 - v[s,j]) <= 0, name=f"battery_charging_ub")
    #                     m.addConstr(B_node[j] - B_node[s] + beta*tij[s,j] - alpha*E_node[s] + C_bat*(1 - v[s,j]) >= 0, name=f"battery_charging_lb")

    # # Battery Capacity constraints
    # # Max Battery
    # for s in F:
    #     m.addConstr(C_bat - B_node[s] - alpha*E_node[s] >= 0, name=f"Max_Battery_Capacity")
    # for i in N:
    #     m.addConstr(B_node[i] - gamma * C_bat >= 0, name=f"Min_Battery_Constraint")
    # m.addConstr(B_node[0] == C_bat, name=f"init_soc[{k}]")

    m.update()
    return m, {"x": x, "v": v, "y": y,  "z": z, "T_node": T_node, "T_veh": T_veh, "D_node": D_node, "B_node": B_node, "E_node": E_node}




