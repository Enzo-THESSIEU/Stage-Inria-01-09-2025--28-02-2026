##### Imports #####

import gurobipy as gb

import Graph_Builder_4 as g

def build_idarp_model_Enhanced_6_Timetabled_Departures(Time_Limit):

    data = g.data   # grab from Graph_builder

    # Unpack sets
    N, P, D, C, R, K = data['N'], data['P'], data['D'], data['C'], data['R'], data['K']
    N_minor, P_minor, D_minor, C_minor = data['N_minor'], data['P_minor'], data['D_minor'], data['C_minor']
    zeroDepot_node, endDepot_node, Departures = data['zeroDepot_node'], data['endDepot_node'], data['Departures']

    Q, T, cij, tij = data['Q'], data['T'], data['cij'], data['tij']
    qr, di, ei, li, ek, lk = data['qr'], data['di'], data['ei'], data['li'], data['ek'], data['lk']
    fi_r, Lbar, pair_pi_di, M, n, n_K = data['fi_r'], data['Lbar'], data['pair_pi_di'], data['M'], data['n'], data['n_K']

    # C_bat, alpha, beta, gamma = data['C_bat'], data['alpha'], data['beta'], data['gamma']

    # Build arc list from tij keys
    A = []
    tij_keys = tij.keys()

    for a in N:
        for b in N:
            if (a[0],b[0]) in tij_keys:
                A.append((a , b))
            if a[0] == b[0]:
                if a[1] != b[1]:
                    A.append((a , b))

    Aset = A
    print("---------------------------------------")
    print("---------------------------------------")
    print("Aset is :", Aset)
    print("---------------------------------------")
    print("---------------------------------------")

    m = gb.Model("IDARP_Enhanced_6_Timetabled_Departures")
    m.setParam('TimeLimit', Time_Limit)

    # Decision variables
    # x_{imjn}^k ∈ {0,1}
    x = m.addVars(K, N, N, vtype=gb.GRB.BINARY, name="x")

    # v_{imjn} ∈ {0,1}
    v = m.addVars(N, N, vtype=gb.GRB.BINARY, name="v")

    # y_{imjn}^r ≥ 0  (implicitly integral in feasible solutions)
    y = m.addVars(R, N, N, vtype=gb.GRB.CONTINUOUS, lb=0, name="y")

    # z_{imjn}^d for (i,j) in C × C 
    z = {}
    for i in C:
        for j in C:
            if i[0] != j[0] : 
                for d in Departures[(i[0],j[0])].keys():
                    z[(d,i,j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")

    # T_im: start of service at node i
    T_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")

    # T_k: time vehicle k leaves the depot
    T_veh  = m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")

    # D_im: time elapsed from vehicle departure until arrival at node i
    D_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

    # # B_im: Battery state at node i.
    # B_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")

    # # E_im: Time spent recharging at station i in F
    # E_node = m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")

    # Objective (1): min sum c_imjn * sum_k x_imjn^k
    m.setObjective(gb.quicksum(cij[(i[0] if isinstance(i, tuple) else i, j[0] if isinstance(j, tuple) else j)] * v[i,j] for (i,j) in Aset), gb.GRB.MINIMIZE)

    # ---------- Constraints ----------

    # (2) Visit pickups/dropoffs exactly once
    for i in P + D:
        m.addConstr(gb.quicksum(v[i, j] for j in N if (i, j) in Aset) == 1, name=f"visit_once_out[{i}]")
        m.addConstr(gb.quicksum(v[j, i] for j in N if (j, i) in Aset) == 1, name=f"visit_once_in[{i}]")

    # (2 annex) Transfer nodes at most once
    for i in C:
        m.addConstr(gb.quicksum(v[i, j] for j in N if (i, j) in Aset) <= 1, name=f"visit_transfer_out[{i}]")
        m.addConstr(gb.quicksum(v[j, i] for j in N if (j, i) in Aset) <= 1, name=f"visit_transfer_in[{i}]")

    # (3) Vehicles leave depot at most once
    for k in K:
        m.addConstr(gb.quicksum(x[k, zeroDepot_node, j] for j in N if (zeroDepot_node, j) in Aset) <= 1,
                    name=f"start_from_depot[{k}]")
    for j in N:
        if (zeroDepot_node, j) in Aset:
            m.addConstr(gb.quicksum(x[k, zeroDepot_node, j] for k in K) == v[zeroDepot_node, j],
                        name=f"link_depot[{j}]")

    # (4) Flow conservation
    for k in K:
        for i in P + D + C:
            m.addConstr(
                gb.quicksum(x[k, j, i] for j in N if (j, i) in Aset) -
                gb.quicksum(x[k, i, j] for j in N if (i, j) in Aset) == 0,
                name=f"flow_x[{k},{i}]"
            )
    for i in P + D + C:
        m.addConstr(
            gb.quicksum(v[j, i] for j in N if (j, i) in Aset) -
            gb.quicksum(v[i, j] for j in N if (i, j) in Aset) == 0,
            name=f"flow_v[{i}]"
        )

    # (5) End at depot
    m.addConstr(gb.quicksum(v[i, endDepot_node] for i in N if (i, endDepot_node) in Aset) == n_K)
    for k in K:
        m.addConstr(gb.quicksum(x[k, i, endDepot_node] for i in N if (i, endDepot_node) in Aset) == 1,
                    name=f"veh_end[{k}]")

    # (6) Passenger balance at transfers
    for r in R:
        for i in C:
            lhs_out = gb.quicksum(y[r, i, j] for j in N if (i, j) in Aset)
            lhs_in  = gb.quicksum(y[r, j, i] for j in N if (j, i) in Aset)
            z_out = gb.quicksum(z[(d, i, j)] for j in C for d in Departures[(i[0], j[0])].keys() if (i, j) in Aset)
            z_in  = gb.quicksum(z[(d, j, i)] for j in C for d in Departures[(j[0], i[0])].keys() if (j, i) in Aset)
            m.addConstr(lhs_out + z_out - lhs_in - z_in == fi_r[(r, i)], name=f"pass_bal_transfer[{r},{i}]")

    # (7) Passenger balance at non-transfers
    for r in R:
        for i in N:
            if i not in C:
                m.addConstr(
                    gb.quicksum(y[r, i, j] for j in N if (i, j) in Aset) -
                    gb.quicksum(y[r, j, i] for j in N if (j, i) in Aset) == fi_r[(r, i)],
                    name=f"pass_bal_node[{r},{i}]"
                )

    # (8) Capacity
    for (i, j) in A:
        m.addConstr(gb.quicksum(qr[r] * y[r, i, j] for r in R) <= Q * v[i, j])
        m.addConstr(gb.quicksum(qr[r] * y[r, i, j] for r in R) <= Q * gb.quicksum(x[k, i, j] for k in K))

    # (9) Timing from depot
    for k in K:
        for j in N:
            if (zeroDepot_node, j) in Aset:
                m.addConstr(T_node[j] >= T_veh[k] + tij[(zeroDepot_node, j)] - M * (1 - x[k, zeroDepot_node, j]))

    # (10) Timing between nodes
    for (i, j) in A:
        m.addConstr(T_node[j] >= T_node[i] + di[i[0]] + tij[(i[0], j[0])] - M * (1 - v[i, j]))
        m.addConstr(T_node[j] >= T_node[i] + di[i[0]] + tij[(i[0], j[0])] - M * (1 - gb.quicksum(x[k, i, j] for k in K)))

    # (11) Fixed-route timing
    for i in C:
        for j in C:
            if (i, j) in Aset:
                for d in Departures[(i[0], j[0])].keys():
                    m.addConstr(T_node[j] >= T_node[i] + tij[(i, j)] - M * (1 - z[(d, i, j)]))

    # (12) Time windows
    for i in P + D:
        m.addConstr(T_node[i] >= ei[i[0]])
        m.addConstr(T_node[i] <= li[i[0]])

    # (13) Ride time
    for i in P:
        drop = pair_pi_di[i]
        if (i, drop) in Aset:
            m.addConstr(T_node[drop] - (T_node[i] + di[i[0]]) >= tij[(i[0], drop[0])])
            m.addConstr(T_node[drop] - (T_node[i] + di[i[0]]) <= Lbar[i[0]])

    # (14) D cumulative time
    for (i, j) in A:
        m.addConstr(D_node[j] >= D_node[i] + di[i[0]] + tij[(i[0], j[0])] - M * (1 - v[i, j]))
        m.addConstr(D_node[j] >= D_node[i] + di[i[0]] + tij[(i[0], j[0])] - M * (1 - gb.quicksum(x[k, i, j] for k in K)))

    # (15) Max duration
    for i in N:
        if (i, endDepot_node) in Aset:
            m.addConstr(D_node[i] + di[i[0]] + tij[(i[0], endDepot_node[0])] <= T)

    # (16) Depot departure window
    for k in K:
        m.addConstr(T_veh[k] >= ek[k])
        m.addConstr(T_veh[k] <= lk[k])

    # (17) Link v and x
    for (i, j) in A:
        if i != zeroDepot_node:
            m.addConstr(v[i, j] == gb.quicksum(x[k, i, j] for k in K))

    # (18) Only one vehicle from depot to j
    for j in P + D + C:
        if (zeroDepot_node, j) in Aset:
            m.addConstr(gb.quicksum(x[k, zeroDepot_node, j] for k in K) == v[zeroDepot_node, j])

    # (19) Subtour elimination: forbid both (i,j) and (j,i)
    for i in N:
        for j in N:
            if i != j and (i, j) in Aset and (j, i) in Aset:
                m.addConstr(v[i, j] + v[j, i] <= 1)
    
    ##### Transfer node strengthening constraints #####
    # (20 - added) A maximum of two transfer nodes can be visited per request
    # for r in R:
    #     m.addConstr(gb.quicksum(v[i,j] for i in C for j in N if (i,j) in Aset) <= 2, name=f"max_transfers[{r}]")

    # (21 - added) If a vehicle visits a transfer node, it must use or have used the fixed-route arc for that request
    for r in R:
        for i in C:
            m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in Aset) == 
                        gb.quicksum(z[(d,i,j)] + z[(d,j,i)] for d in Departures[i[0],j[0]].keys() for j in C if (i,j) in z),
                        name=f"transfer_sync[{r},{i}]")
            
    ##### Scheduled PT Constraints #####
    for i in N:
        if i not in C:
            m.addConstr(M*(1 - gb.quicksum(gb.quicksum(z[(d,i,j)] for d in Departures[i[0],j[0]] for j in C if (i,j) in Aset))) + 
                        gb.quicksum(gb.quicksum(Departures[i[0],j[0]][d] * z[(d,i,j)] for d in Departures[i[0],j[0]].keys() for j in C if (i,j) in Aset)) - T_node[i] >= 0 , name="PT Departure Time Constraint")

    for i in C:
        m.addConstr(T_node[i] - gb.quicksum(gb.quicksum((Departures[i[0],j[0]][d] + tij[(i[0],j[0])]) * z[(d,i,j)] for d in Departures[i[0],j[0]].keys() for j in C if (i,j) in Aset)) >= 0 , name="PT Arrival Time Constraint")


    m.update()
    return m, {"x": x, "v": v, "y": y,  "z": z, "T_node": T_node, "T_veh": T_veh, "D_node": D_node} # "B_node": B_node, "E_node": E_node}




