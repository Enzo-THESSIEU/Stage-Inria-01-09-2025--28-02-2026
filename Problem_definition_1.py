##### Imports #####

import gurobipy as gb
import Graph_builder_1 as g

##### Original Function Definition #####

def build_idarp_model_Original(Time_Limit):

    data = g.data_original   # grab from Graph_builder

    # Unpack sets
    N, P, D, C, R, K, Cr = data['N'], data['P'], data['D'], data['C'], data['R'], data['K'], data['Cr']
    Q, T, cij, tij, qr, di, ei, li, ek, lk = data['Q'], data['T'], data['cij_original'], data['tij_original'], data['qr'], data['di'], data['ei'], data['li'], data['ek'], data['lk']
    fi_r, Lbar, endDepot, zeroDepot, pair_pi_di, M = data['fi_r'], data['Lbar'], data['endDepot'], data['zeroDepot'], data['pair_pi_di'], data['M']

    m = gb.Model("IDARP_Original")
    m.setParam('TimeLimit', Time_Limit) 

    # Decision variables
    # x_{ij}^k ∈ {0,1}
    x = m.addVars(K, N, N, vtype=gb.GRB.BINARY, name="x")

    # y_{ij}^r ≥ 0  (implicitly integral in feasible solutions)
    y = m.addVars(R, N, N, vtype=gb.GRB.CONTINUOUS, lb=0.0, name="y")

    # z_{ij} for (i,j) in C_r × C_r — we index it as z[r,i,j] to avoid ambiguity
    z = {}
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                z[(i,j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{r},{i},{j}]")

    # B_i: start of service at node i
    B_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")

    # B_k: time vehicle k leaves the depot
    B_veh  = m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="B_veh")

    # D_i: time elapsed from vehicle departure until arrival at node i
    D_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

    # Objective (1): min sum c_ij * sum_k x_ij^k
    m.setObjective(gb.quicksum(cij[i,j] * x[k,i,j] for k in K for i in N for j in N if (i,j) in cij), gb.GRB.MINIMIZE)

    # ---------- Constraints ----------

    # (2) Every pickup/dropoff node visited exactly once by some vehicle
    for i in P + D:
        m.addConstr(gb.quicksum(x[k,i,j] for k in K for j in N if (i,j) in tij) == 1, name=f"visit_once[{i}]")

    # (3) Each vehicle leaves the origin depot exactly once
    for k in K:
        m.addConstr(gb.quicksum(x[k,zeroDepot,j] for j in N if (zeroDepot,j) in tij) == 1, name=f"start_from_depot[{k}]")

    # (4) Vehicle flow conservation at pickups, dropoffs, and transfers
    for k in K:
        for i in (P + D + C):
            m.addConstr(gb.quicksum(x[k,j,i] for j in N if (j,i) in tij) - gb.quicksum(x[k,i,j] for j in N if (i,j) in tij) == 0,
                        name=f"flow_cons[{k},{i}]")

    # (5) Each vehicle ends at the destination depot exactly once
    for k in K:
        m.addConstr(gb.quicksum(x[k,i,endDepot] for i in N if (i,endDepot) in tij) == 1, name=f"end_at_depot[{k}]")

    # (6) Passenger balance at transfer nodes (can be via DRT or fixed route)
    for r in R:
        for i in Cr[r]:
            m.addConstr(
                gb.quicksum(y[r,i,j] for j in N if (i,j) in tij) + gb.quicksum(z[(i,j)] for j in Cr[r] if (i,j) in tij)
                - gb.quicksum(y[r,j,i] for j in N if (j,i) in tij) - gb.quicksum(z[(j,i)] for j in Cr[r]if (j,i) in tij)
                == fi_r[(r,i)],
                name=f"pass_bal_transfer[{r},{i}]"
            )

    # (7) Passenger balance at non-transfer nodes (only DRT arcs)
    for r in R:
        for i in N:
            if i not in Cr[r]:
                m.addConstr(
                    gb.quicksum(y[r,i,j] for j in N if (i,j) in tij) - gb.quicksum(y[r,j,i] for j in N if (j,i) in tij) == fi_r[(r,i)],
                    name=f"pass_bal_nontransfer[{r},{i}]"
                )

    # (8) Capacity on vehicle arcs
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * gb.quicksum(x[k,i,j] for k in K),
                            name=f"capacity[{i},{j}]")

    # (9) Timing from depot to first served node on vehicle k
    for k in K:
        for j in N:
            if (i,j) in tij and (zeroDepot,j) in tij:
                m.addConstr(
                    B_node[j] >= B_veh[k] + tij[zeroDepot,j] - M*(1 - x[k,zeroDepot,j]),
                    name=f"time_dep_first[{k},{j}]"
                )

    # (10) Timing between consecutive service nodes (DRT)
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(
                    B_node[j] >= B_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                    name=f"time_link_drt[{i},{j}]"
                )

    # (11) Timing along fixed-route between transfer nodes for request r
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                if (i,j) in tij:
                    m.addConstr(
                        B_node[j] >= B_node[i] + tij[i,j] - M*(1 - z[(i,j)]),
                        name=f"time_link_fixed[{r},{i},{j}]"
                    )

    # (12) Time windows at pickup and dropoff nodes
    ### for i in (P + D): ### No time windows for dropoffs in this example
    for i in P:
        m.addConstr(B_node[i] >= ei[i], name=f"tw_lo[{i}]")
        m.addConstr(B_node[i] <= li[i], name=f"tw_hi[{i}]")

    # (13) Ride-time constraints: t_{i,n+i} <= B_{n+i} - (B_i + d_i) <= Lbar_i  for pickups i
    for i in P:
        drop = pair_pi_di[i]  # dropoff node for pickup i (usually n+i)
        if (i,drop) in tij:
            m.addConstr(B_node[drop] - (B_node[i] + di[i]) >= tij[i,drop], name=f"ride_lo[{i}]")
            m.addConstr(B_node[drop] - (B_node[i] + di[i]) <= Lbar[i],   name=f"ride_hi[{i}]")

    # (14) Cumulative time variable D (mirrors (10) without B, used for duration in (15))
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(
                    D_node[j] >= D_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                    name=f"D_link[{i},{j}]"
                )

    # (15) Vehicle maximum route duration T (ensure from node i to end depot fits in T)
    for i in N:
        if (i,endDepot) in tij:
            m.addConstr(D_node[i] + di[i] + tij[i,endDepot] <= T, name=f"duration[{i}]")

    # (16) Depot departure window for each vehicle
    for k in K:
        m.addConstr(B_veh[k] >= ek[k], name=f"veh_dep_lo[{k}]")
        m.addConstr(B_veh[k] <= lk[k], name=f"veh_dep_hi[{k}]")

    # (17) x are binary    — already set in var types
    # (18) y >= 0          — already set in var types
    # (19) z are binary    — already set in var types

    m.update()
    return m, {"x": x, "y": y, "z": z, "B_node": B_node, "B_veh": B_veh, "D_node": D_node}

##### First Enhancement with Arc Elimination #####

def build_idarp_model_Enhanced_1(Time_Limit):

    data = g.data   # grab from Graph_builder

    # Unpack sets
    N, P, D, C, R, K, Cr = data['N'], data['P'], data['D'], data['C'], data['R'], data['K'], data['Cr']
    Q, T, cij, tij, qr, di, ei, li, ek, lk = data['Q'], data['T'], data['cij'], data['tij'], data['qr'], data['di'], data['ei'], data['li'], data['ek'], data['lk']
    fi_r, Lbar, endDepot, zeroDepot, pair_pi_di, M = data['fi_r'], data['Lbar'], data['endDepot'], data['zeroDepot'], data['pair_pi_di'], data['M']

    m = gb.Model("IDARP_Enhanced_1")
    m.setParam('TimeLimit', Time_Limit)

    # Decision variables
    # x_{ij}^k ∈ {0,1}
    x = m.addVars(K, N, N, vtype=gb.GRB.BINARY, name="x")

    # y_{ij}^r ≥ 0  (implicitly integral in feasible solutions)
    y = m.addVars(R, N, N, vtype=gb.GRB.CONTINUOUS, lb=0.0, name="y")

    # z_{ij} for (i,j) in C_r × C_r — we index it as z[r,i,j] to avoid ambiguity
    z = {}
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                z[(i,j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{r},{i},{j}]")

    # B_i: start of service at node i
    B_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")

    # B_k: time vehicle k leaves the depot
    B_veh  = m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="B_veh")

    # D_i: time elapsed from vehicle departure until arrival at node i
    D_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

    # Objective (1): min sum c_ij * sum_k x_ij^k
    m.setObjective(gb.quicksum(cij[i,j] * x[k,i,j] for k in K for i in N for j in N if (i,j) in cij), gb.GRB.MINIMIZE)

    # ---------- Constraints ----------

    # (2) Every pickup/dropoff node visited exactly once by some vehicle
    for i in P + D:
        m.addConstr(gb.quicksum(x[k,i,j] for k in K for j in N if (i,j) in tij) == 1, name=f"visit_once[{i}]")

    # (3) Each vehicle leaves the origin depot exactly once
    for k in K:
        m.addConstr(gb.quicksum(x[k,zeroDepot,j] for j in N if (zeroDepot,j) in tij) == 1, name=f"start_from_depot[{k}]")

    # (4) Vehicle flow conservation at pickups, dropoffs, and transfers
    for k in K:
        for i in (P + D + C):
            m.addConstr(gb.quicksum(x[k,j,i] for j in N if (j,i) in tij) - gb.quicksum(x[k,i,j] for j in N if (i,j) in tij) == 0,
                        name=f"flow_cons[{k},{i}]")

    # (5) Each vehicle ends at the destination depot exactly once
    for k in K:
        m.addConstr(gb.quicksum(x[k,i,endDepot] for i in N if (i,endDepot) in tij) == 1, name=f"end_at_depot[{k}]")

    # (6) Passenger balance at transfer nodes (can be via DRT or fixed route)
    for r in R:
        for i in Cr[r]:
            m.addConstr(
                gb.quicksum(y[r,i,j] for j in N if (i,j) in tij) + gb.quicksum(z[(i,j)] for j in Cr[r] if (i,j) in tij)
                - gb.quicksum(y[r,j,i] for j in N if (j,i) in tij) - gb.quicksum(z[(j,i)] for j in Cr[r]if (j,i) in tij)
                == fi_r[(r,i)],
                name=f"pass_bal_transfer[{r},{i}]"
            )

    # (7) Passenger balance at non-transfer nodes (only DRT arcs)
    for r in R:
        for i in N:
            if i not in Cr[r]:
                m.addConstr(
                    gb.quicksum(y[r,i,j] for j in N if (i,j) in tij) - gb.quicksum(y[r,j,i] for j in N if (j,i) in tij) == fi_r[(r,i)],
                    name=f"pass_bal_nontransfer[{r},{i}]"
                )

    # (8) Capacity on vehicle arcs
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * gb.quicksum(x[k,i,j] for k in K),
                            name=f"capacity[{i},{j}]")

    # (9) Timing from depot to first served node on vehicle k
    for k in K:
        for j in N:
            if (i,j) in tij and (zeroDepot,j) in tij:
                m.addConstr(
                    B_node[j] >= B_veh[k] + tij[zeroDepot,j] - M*(1 - x[k,zeroDepot,j]),
                    name=f"time_dep_first[{k},{j}]"
                )

    # (10) Timing between consecutive service nodes (DRT)
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(
                    B_node[j] >= B_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                    name=f"time_link_drt[{i},{j}]"
                )

    # (11) Timing along fixed-route between transfer nodes for request r
    for r in R:
        for i in Cr[r]:
            for j in Cr[r]:
                if (i,j) in tij:
                    m.addConstr(
                        B_node[j] >= B_node[i] + tij[i,j] - M*(1 - z[(i,j)]),
                        name=f"time_link_fixed[{r},{i},{j}]"
                    )

    # (12) Time windows at pickup and dropoff nodes
    ### for i in (P + D): ### No time windows for dropoffs in this example
    for i in P + D:
        m.addConstr(B_node[i] >= ei[i], name=f"tw_lo[{i}]")
        m.addConstr(B_node[i] <= li[i], name=f"tw_hi[{i}]")

    # (13) Ride-time constraints: t_{i,n+i} <= B_{n+i} - (B_i + d_i) <= Lbar_i  for pickups i
    for i in P:
        drop = pair_pi_di[i]  # dropoff node for pickup i (usually n+i)
        if (i,drop) in tij:
            m.addConstr(B_node[drop] - (B_node[i] + di[i]) >= tij[i,drop], name=f"ride_lo[{i}]")
            m.addConstr(B_node[drop] - (B_node[i] + di[i]) <= Lbar[i],   name=f"ride_hi[{i}]")

    # (14) Cumulative time variable D (mirrors (10) without B, used for duration in (15))
    for i in N:
        for j in N:
            if (i,j) in tij:
                m.addConstr(
                    D_node[j] >= D_node[i] + di[i] + tij[i,j] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
                    name=f"D_link[{i},{j}]"
                )

    # (15) Vehicle maximum route duration T (ensure from node i to end depot fits in T)
    for i in N:
        if (i,endDepot) in tij:
            m.addConstr(D_node[i] + di[i] + tij[i,endDepot] <= T, name=f"duration[{i}]")

    # (16) Depot departure window for each vehicle
    for k in K:
        m.addConstr(B_veh[k] >= ek[k], name=f"veh_dep_lo[{k}]")
        m.addConstr(B_veh[k] <= lk[k], name=f"veh_dep_hi[{k}]")

    # (17) x are binary    — already set in var types
    # (18) y >= 0          — already set in var types
    # (19) z are binary    — already set in var types

    m.update()
    return m, {"x": x, "y": y, "z": z, "B_node": B_node, "B_veh": B_veh, "D_node": D_node}

# m, vars_ = build_idarp_model_Original(5)
# m.optimize()
# print("Obj:", m.ObjVal)

# # x0 = vars_["x"]
# # print({ (k,i,j): var.X for (k,i,j), var in x0.items() })

# # Problem_definition_1.py (after m.optimize())
# # if m.status == gb.GRB.OPTIMAL:
# print("Obj:", m.ObjVal)


# used_arcs_vehicle_1, used_arcs_vehicle_2 = extract_routes(
#     vars_,
#     nodes=g.data['nodes'],
#     N=g.data['N'],
#     K=g.data['K'],
#     zeroDepot=g.data['zeroDepot'],
#     endDepot=g.data['endDepot']
# )

# print("\nVehicle 1 route:", used_arcs_vehicle_1)
# print("Vehicle 2 route:", used_arcs_vehicle_2)

# print(g.data['fi_r'][3,7])
# print(g.data['fi_r'][1,7])
# print(g.data['fi_r'][4,4])

# r, i = 3, 7  # request 3, node 7
# val_1 = sum(vars_['y'][r, i, j].X for j in g.data['N'] if (i, j) in g.data['tij'])
# print(f"Sum of y[{r},{i},j] over all j = {val_1}")
# val_2 = sum(vars_['y'][r, j, i].X for j in g.data['N'] if (j, i) in g.data['tij'])
# print(f"Sum of y[{r},j,{i}] over all j = {val_2}")

# # else:
# #     print(f"Optimization ended with status: {m.status}")


