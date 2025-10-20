# ########## Parametres.py ##########


# # ##### Necessary Constants #####
# # max_visits_transfer = 10

# # ##### Build Nodes/Arcs and Parametres #####

# # ### === FUNCTIONS ===

# # def base(x):
# #     return x[0] if isinstance(x, tuple) else x


# # def build_t_transfer(t, n_requests, n_vehicles, n_trans_nodes, duplicate_transfers=True, ev_constraints=False):
# #     """
# #     Build the extended travel-time matrix including transfer and/or charging nodes.

# #     Parameters
# #     ----------
# #     t : np.ndarray
# #         Base travel time matrix (n_base × n_base).
# #     n_requests : int
# #         Number of requests.
# #     n_vehicles : int
# #         Number of vehicles.
# #     n_trans_nodes : int
# #         Number of transfer nodes per cluster.
# #     duplicate_transfers : bool, optional
# #         If True, each request gets its own duplicated set of transfer nodes.
# #         If False, transfers are shared across requests.
# #     ev_constraints : bool, optional
# #         If True, also add charging station nodes for each vehicle.

# #     Returns
# #     -------
# #     t_transfer : np.ndarray
# #         The full expanded travel-time matrix including transfers and (optionally) charging stations.
# #     """

# #     n_base = t.shape[0]                # Base node count (P + D + depots, etc.)
# #     n_transfers = n_trans_nodes * (n_requests if duplicate_transfers else 1)
# #     n_charging = n_vehicles * n_trans_nodes if ev_constraints else 0
# #     n_total = n_base + n_transfers + n_charging

# #     t_transfer = np.zeros((n_total, n_total))

# #     # ---- Base block ----
# #     t_transfer[:n_base, :n_base] = t

# #     # Index ranges for clarity
# #     trans_start = n_base
# #     trans_end = trans_start + n_transfers
# #     charge_start = trans_end
# #     charge_end = charge_start + n_charging

# #     # ---- Base → Transfers ----
# #     if n_transfers > 0:
# #         base_to_trans = np.tile(t[:, 10:10 + n_trans_nodes],
# #                         (1, (n_requests if duplicate_transfers else 1)))

# #         trans_to_base = np.tile(t[10:10 + n_trans_nodes, :],
# #                                 ((n_requests if duplicate_transfers else 1), 1))

# #         trans_block = np.tile(t[10:10 + n_trans_nodes, 10:10 + n_trans_nodes],
# #                             ((n_requests if duplicate_transfers else 1),
# #                             (n_requests if duplicate_transfers else 1)))
# #         t_transfer[trans_start:trans_end, trans_start:trans_end] = trans_block

# #     # ---- Base ↔ Charging stations ----
# #     if ev_constraints and n_charging > 0:
# #         base_to_charge = np.tile(t[:, 10:10 + n_trans_nodes], (1, n_vehicles))
# #         charge_to_base = np.tile(t[10:10 + n_trans_nodes, :], (n_vehicles, 1))
# #         charge_block = np.tile(t[10:10 + n_trans_nodes, 10:10 + n_trans_nodes],
# #                                (n_vehicles, n_vehicles))

# #         # Base → Chargers
# #         t_transfer[:n_base, charge_start:charge_end] = base_to_charge
# #         # Chargers → Base
# #         t_transfer[charge_start:charge_end, :n_base] = charge_to_base
# #         # Chargers ↔ Chargers
# #         t_transfer[charge_start:charge_end, charge_start:charge_end] = charge_block

# #         # Optionally connect transfers ↔ chargers (only if both exist)
# #         if n_transfers > 0:
# #             t_transfer[trans_start:trans_end, charge_start:charge_end] = np.mean(charge_block)
# #             t_transfer[charge_start:charge_end, trans_start:trans_end] = np.mean(charge_block)

# #     return t_transfer

# # def build_nodes(base_nodes, n_requests, n_vehicles, n_trans_nodes, duplicate_nodes=True, ev_constraints=True):
# #     """
# #     Build node list with consistent numbering:
# #       - Transfer nodes start at 10, 11, 12, ...
# #       - Charging station nodes start right after the last transfer node.

# #     Parameters
# #     ----------
# #     base_nodes : np.ndarray
# #         Array of base nodes (depot, pickups, dropoffs)
# #     n_requests : int
# #         Number of requests
# #     n_vehicles : int
# #         Number of vehicles
# #     n_trans_nodes : int
# #         Number of physical transfer nodes per request or cluster
# #     duplicate_nodes : bool
# #         Whether to duplicate transfers for each request
# #     ev_constraints : bool
# #         Whether to include charging station nodes

# #     Returns
# #     -------
# #     nodes : np.ndarray
# #         Full node array including base, transfer, and optional charging nodes
# #     """
# #     transfers, charging_stations = [], []

# #     # --- TRANSFER NODES ---
# #     transfer_start_id = 10
# #     if duplicate_nodes:
# #         # One set of transfers per request
# #         for r in range(1, n_requests + 1):
# #             for t_id in range(n_trans_nodes):
# #                 node_id = transfer_start_id + (r - 1) * n_trans_nodes + t_id
# #                 transfers.append([node_id, "transfer", str(r), "", "", 5, 0])
# #     else:
# #         # Single shared set of transfers
# #         for t_id in range(n_trans_nodes):
# #             node_id = transfer_start_id + t_id
# #             transfers.append([node_id, "transfer", "", "", "", 5, 0])

# #     # --- CHARGING STATIONS ---
# #     if ev_constraints:
# #         charge_start_id = transfer_start_id + n_trans_nodes * (n_requests if duplicate_nodes else 1)
# #         if duplicate_nodes:
# #             for k in range(1, n_vehicles + 1):
# #                 for c_id in range(n_trans_nodes):
# #                     node_id = charge_start_id + (k - 1) * n_trans_nodes + c_id
# #                     charging_stations.append([node_id, "charging station", "", 0, 86400, 0, 0])
# #         else:
# #             for c_id in range(n_trans_nodes):
# #                 node_id = charge_start_id + c_id
# #                 charging_stations.append([node_id, "charging station", "", 0, 86400, 0, 0])

# #     # --- SAFE CONVERSION ---
# #     if len(transfers) == 0:
# #         transfers = np.empty((0, base_nodes.shape[1]), dtype=object)
# #     if len(charging_stations) == 0:
# #         charging_stations = np.empty((0, base_nodes.shape[1]), dtype=object)

# #     # --- FINAL STACK ---
# #     return np.vstack([
# #         base_nodes,
# #         np.array(transfers, dtype=object),
# #         np.array(charging_stations, dtype=object)
# #     ])

# # def build_node_id_sets(nodes, ev_constraints=False):
# #     N = nodes[:, 0].astype(int).tolist()
# #     P = nodes[nodes[:, 1] == "pickup", 0].astype(int).tolist()
# #     D = nodes[nodes[:, 1] == "dropoff", 0].astype(int).tolist()
# #     C = nodes[nodes[:, 1] == "transfer", 0].astype(int).tolist()
# #     zeroDepot = int(nodes[nodes[:, 1] == "depot0", 0][0])
# #     endDepot = int(nodes[nodes[:, 1] == "depot1", 0][0])
# #     F = []
# #     if ev_constraints:
# #         F = nodes[nodes[:, 1] == "charging station", 0].astype(int).tolist()
# #     return N, P, D, C, F, zeroDepot, endDepot

# # def build_node_imjn_sets(nodes, max_visits_transfer, ev_constraints=False):
# #     N_minor, P_minor, D_minor, C_minor, F_minor, zeroDepot_node, endDepot_node = build_node_id_sets(nodes)
# #     # Wrap depots as (id,1)
# #     zeroDepot_node = (zeroDepot_node, 1)
# #     endDepot_node = (endDepot_node, 1)

# #     # Expand pickups/dropoffs as (i,1)
# #     P = [(i, 1) for i in P_minor]
# #     D = [(i, 1) for i in D_minor]

# #     # Expand transfer nodes with artificial visits
# #     C = [(i, m) for m in range(max_visits_transfer) for i in C_minor]
# #     F = []
# #     if ev_constraints:
# #         F = [(i, m) for i in (F_minor or []) for m in range(max_visits_transfer)]

# #     # All nodes in tuple format
# #     N = [zeroDepot_node, endDepot_node] + P + D + C + F

# #     return N, P, D, C, F, zeroDepot_node, endDepot_node

# # def build_requests(nodes, use_imjn=False):
# #     if use_imjn:
# #         r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "transfer"]), 2]
# #         r_values = r_values[r_values != ""]
# #         R = sorted(set(r_values.astype(int)))
# #         return R
# #     else:
# #         r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "transfer"]), 2]
# #         r_values = r_values[r_values != ""]
# #         R = sorted(set(r_values.astype(int)))
# #         Cr = {r: nodes[(nodes[:, 1] == "transfer") & (nodes[:, 2] == str(r)), 0].astype(int).tolist()
# #             for r in R}
# #         return R, Cr

# # def build_service_params(nodes, t, P, D):
# #     di = {int(row[0]): float(row[5]) for row in nodes if row[5] != ""}
# #     ei = {int(row[0]): float(row[3]) for row in nodes if row[3] != ""}
# #     li = {int(row[0]): float(row[4]) for row in nodes if row[4] != ""}
# #     r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff"]), 2]
# #     r_values = r_values[r_values != ""]
# #     R = sorted(set(r_values.astype(int)))
# #     qr = {int(r): 1 for r in R}
# #     tij = {(i, j): float(t[i, j]) for i in range(len(t)) for j in range(len(t))}
# #     pair_pi_di = {p: d for p, d in zip(P, D)}
# #     Lbar = {self.base(i): 1800 for i in P}
# #     for p in pair_pi_di.keys():
# #         d = pair_pi_di[p]
# #         pid = p[0] if isinstance(p, tuple) else p
# #         did = d[0] if isinstance(d, tuple) else d
# #         e = ei[pid] + di[pid] + tij[pid, did]
# #         l = li[pid] + di[pid] + tij[pid, did]
# #         ei[did] = e
# #         li[did] = l
# #     return di, ei, li, qr, pair_pi_di, Lbar

# # def build_fi_r(nodes, N, R, P, D):
# #     fi_r = {(int(r), i): 0 for r in R for i in N}
# #     for p, d in zip(P, D):
# #         r = int(nodes[nodes[:, 0] == self.base(p), 2][0])
# #         fi_r[r, p] = 1
# #         fi_r[r, d] = -1
# #     return fi_r

# # def build_departures(t, C_minor, interval=5, planning_horizon=24*60):
# #     """
# #     Build Departures dictionary with forward and backward arcs (including skip arcs).

# #     Args:
# #         t (np.array): travel time matrix (base or expanded).
# #         C_minor (list): list of transfer node ids (integers).
# #         interval (int): departure interval in seconds (default 300 = 5 minutes).
# #         planning_horizon (int): total time horizon (default 24h).

# #     Returns:
# #         dict: Departures[(i,j)] = {a: departure_time}
# #     """
# #     Departures = {}
# #     n_intervals = planning_horizon // interval
# #     transfer_nodes = sorted(C_minor)

# #     # === Forward direction (left → right) ===
# #     left_terminal = transfer_nodes[0]
# #     for j in range(1, len(transfer_nodes)):
# #         Departures[(left_terminal, transfer_nodes[j])] = {
# #             a: interval * a for a in range(int(n_intervals))
# #         }

# #     # Forward arcs (including skips)
# #     for i in range(1, len(transfer_nodes)):
# #         for j in range(i+1, len(transfer_nodes)):
# #             Departures[(transfer_nodes[i], transfer_nodes[j])] = {
# #                 a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k+1]] for k in range(i, j)))
# #                 for a, d in Departures[(left_terminal, transfer_nodes[i])].items()
# #             }

# #     # === Backward direction (right → left) ===
# #     right_terminal = transfer_nodes[-1]
# #     for j in range(len(transfer_nodes)-1):
# #         Departures[(right_terminal, transfer_nodes[j])] = {
# #             a: interval * a for a in range(int(n_intervals))
# #         }

# #     # Backward arcs (including skips)
# #     for i in range(len(transfer_nodes)-2, -1, -1):
# #         for j in range(i-1, -1, -1):
# #             Departures[(transfer_nodes[i], transfer_nodes[j])] = {
# #                 a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k-1]] for k in range(i, j, -1)))
# #                 for a, d in Departures[(right_terminal, transfer_nodes[i])].items()
# #             }

# #     # Departures[(10,10)] = {a + 1:300*a for a in range(int(n_intervals))}
# #     # Departures[(12,12)] = {a + 1:300*a for a in range(int(n_intervals))}

# #     return Departures

# # ##### Arc Elimination #####

# # def build_tij(t, nodes, n_requests, n_vehicles, n_trans_nodes,
# #               include_charging=True, duplicate=True, eliminate_arcs=False,
# #               di=None, ei=None, li=None, Lbar=None, P=None, D=None, N=None,
# #               zeroDepot=None, endDepot=None, Cr=None, n=None, pair_pi_di=None):
    
# #     """
# #     Build tij and cij with options for charging, duplication, and arc elimination.

# #     Args:
# #         t (np.array): base travel time matrix
# #         nodes (np.array): nodes array
# #         n_requests (int): number of requests
# #         n_vehicles (int): number of vehicles
# #         n_trans_nodes (int): number of physical transfer nodes
# #         include_charging (bool): include charging stations
# #         duplicate (bool): duplicate transfer and charging nodes per request/vehicle
# #         eliminate_arcs (bool): apply arc elimination rules
# #         di, ei, li, Lbar, P, D, N, zeroDepot, endDepot, Cr, n: optional sets/params 
# #             needed for arc elimination.

# #     Returns:
# #         tij (dict): travel time dictionary
# #         cij (dict): cost dictionary
# #     """

# #     # === Step 1. Build expanded t_transfer ===
# #     n_base = t.shape[0]
# #     n_trans = n_trans_nodes * (n_requests - 1 if duplicate else 1)
# #     n_charge = n_vehicles * n_trans_nodes if include_charging else (n_trans_nodes if not duplicate else 0)
# #     n_total = n_base + n_trans + n_charge

# #     t_transfer = np.zeros((n_total, n_total))

# #     # Base block
# #     t_transfer[:n_base, :n_base] = t

# #     if n_trans > 0 or n_charge > 0:
# #         factor = (n_trans // n_trans_nodes + n_charge // n_trans_nodes)
# #         # base → transfer/charging
# #         t_transfer[:n_base, n_base:] = np.tile(t[:, 10:10+n_trans_nodes], (1, factor))
# #         # transfer/charging → base
# #         t_transfer[n_base:, :n_base] = np.tile(t[10:10+n_trans_nodes, :], (factor, 1))
# #         # transfer/charging → transfer/charging
# #         t_transfer[n_base:, n_base:] = np.tile(t[10:10+n_trans_nodes, 10:10+n_trans_nodes], (factor, factor))

# #     # === Step 2. Build tij, cij dicts ===
# #     tij = {(self.base(i), self.base(j)): float(t_transfer[self.base(i), self.base(j)]) for i in N for j in N}

# #     # === Step 3. Arc elimination (optional) ===
# #     if eliminate_arcs and all(v is not None for v in [di, ei, li, Lbar, P, D, N, zeroDepot, endDepot, n]):

# #         # Depot restrictions
# #         for i in N:
# #             if (i, zeroDepot) in tij:
# #                 del tij[(self.base(i), self.base(zeroDepot))]
# #             if (endDepot, i) in tij:
# #                 del tij[(self.base(endDepot), self.base(i))]

# #         for i in D:
# #             if (zeroDepot, i) in tij:
# #                 del tij[(self.base(zeroDepot), self.base(i))]
# #         for i in P:
# #             if (i, endDepot) in tij:
# #                 del tij[(self.base(i), self.base(endDepot))]

# #         # Remove loops
# #         for i in N:
# #             if (i, i) in tij:
# #                 del tij[(self.base(i), self.base(i))]

# #         # Drop-off → pickup arcs
# #         for p, d in pair_pi_di.items():
# #             if (d, p) in tij:
# #                 del tij[(self.base(d), self.base(p))]

# #         # Time-window infeasibility
# #         for i in N:
# #             for j in N:
# #                 if (i, j) in tij and i in ei and i in di and j in li:
# #                     if ei[self.base(i)] + di[self.base(i)] + tij[(self.base(i), self.base(j))] >= li[self.base(j)]:
# #                         del tij[(self.base(i), self.base(j))]

# #         # Ride time infeasibility
# #         for i in P:
# #             for j in N:
# #                 if (self.base(i), self.base(j)) in tij and (self.base(j), self.base(pair_pi_di[i])) in tij:
# #                     if tij[(self.base(i), self.base(j))] + di[self.base(j)] + tij[(self.base(j), self.base(pair_pi_di[i]))] >= Lbar[self.base(i)]:
# #                         del tij[(self.base(i), self.base(j))]
# #                         del tij[(self.base(j), self.base(pair_pi_di[i]))]

# #         # Path infeasibility
# #         for i in P:
# #             for j in P:
# #                 if i != j:
# #                     if (self.base(j), self.base(i)) in tij and (self.base(i), self.base(pair_pi_di[j])) in tij and (self.base(pair_pi_di[j]), self.base(pair_pi_di[i])) in tij:
# #                         if tij[(self.base(j), self.base(i))] + di[self.base(i)] + tij[(self.base(i), self.base(pair_pi_di[i]))] + di[self.base(pair_pi_di[i])] + tij[(self.base(pair_pi_di[j]), self.base(pair_pi_di[i]))] > Lbar[self.base(i)]:
# #                             del tij[(self.base(i), self.base(pair_pi_di[i]))]
# #         for i in P:
# #             for j in P:
# #                 if i != j:
# #                     if (self.base(i), self.base(pair_pi_di[i])) in tij and (self.base(pair_pi_di[i]), self.base(j)) in tij and (self.base(j), self.base(pair_pi_di[i])) in tij:
# #                         if tij[(self.base(i), self.base(pair_pi_di[i]))] + di[self.base(pair_pi_di[i])] + tij[(self.base(pair_pi_di[i]), self.base(j))] + di[self.base(j)] + tij[(self.base(j), self.base(pair_pi_di[i]))] > Lbar[self.base(i)]:
# #                             del tij[(self.base(pair_pi_di[i]), self.base(j))]

# #         # Transfer-specific cleaning only if Cr exists
# #         if Cr is not None:
# #             for i in P:
# #                 drop_node = pair_pi_di[i]
# #                 if i in Cr:
# #                     for Si in Cr[i]:
# #                         if (drop_node, Si) in tij:
# #                             del tij[(self.base(drop_node), self.base(Si))]
# #                         if (Si, i) in tij:
# #                             del tij[(self.base(Si), self.base(i))]

# #     cij = tij.copy()
# #     return tij, cij

# # ### === USAGE ===

# # def pack_data(duplicate_transfers = True, arc_elimination = True, ev_constraints = False, use_imjn = False):

# #     # Parametres

# #     ### Base travel time matrix
# #     t = np.array([
# #         [  0,  23, 173, 131, 195, 170, 129, 160, 153,   0, 120,  21, 163],
# #         [ 23,   0, 166, 123, 173, 179, 125, 138, 168,  23, 135,  26, 159],
# #         [173, 166,   0,  57, 289,  41, 268, 254,  79, 173,  57, 191, 301],
# #         [131, 123,  57,   0, 258,  72, 236, 222,  60, 131,  29, 148, 270],
# #         [195, 173, 289, 258,   0, 318,  73,  35, 306, 195, 273, 174,  39],
# #         [170, 179,  41,  72, 318,   0, 296, 282,  53, 170,  51, 191, 329],
# #         [129, 125, 268, 236,  73, 296,   0,  48, 282, 129, 248, 108,  34],
# #         [160, 138, 254, 222,  35, 282,  48,   0, 271, 160, 237, 139,  57],
# #         [153, 168,  79,  60, 306,  53, 282, 271,   0, 153,  35, 174, 316],
# #         [  0,  23, 173, 131, 195, 170, 129, 160, 153,   0, 120,  21, 163],
# #         [120, 135,  57,  29, 273,  51, 248, 237,  35, 120,   0, 141, 283],
# #         [ 21,  26, 191, 148, 174, 191, 108, 139, 174,  21, 141,   0, 142],
# #         [163, 159, 301, 270,  39, 329,  34,  57, 316, 163, 283, 142,   0]
# #     ], dtype=float)

# #     ### Base nodes
# #     nodes = np.array([
# #         [0, "depot0", "", 0, 86400, 0, 0],
# #         [9, "depot1", "", 0, 86400, 0, 0],
# #         [1, "pickup", 1, 50, 950, 5, 1],
# #         [5, "dropoff", 1, "", "", 5, 0],
# #         [2, "pickup", 2, 350, 1250, 5, 1],
# #         [6, "dropoff", 2, "", "", 5, 0],
# #         [3, "pickup", 3, 650, 1550, 5, 1],
# #         [7, "dropoff", 3, "", "", 5, 0],
# #         [4, "pickup", 4, 950, 1850, 5, 1],
# #         [8, "dropoff", 4, "", "", 5, 0],
# #     ], dtype=object)

# #     n_requests = 4
# #     n_vehicles = 2 
# #     n_trans_nodes = 3

# #     Departures = None
# #     Cr = None

# #     # Build matrices and nodes
# #     t_transfer = build_t_transfer(t, n_requests, n_vehicles, n_trans_nodes, duplicate_transfers=duplicate_transfers, ev_constraints=ev_constraints)
# #     # print("size of t_transfer: ", np.size(t_transfer))
# #     nodes = build_nodes(nodes, n_requests, n_vehicles, n_trans_nodes, self.duplicate_nodes=duplicate_transfers, self.ev_constraints=ev_constraints)
# #     # print(nodes)

# #     # Build sets Id
# #     # N, P, D, C, F, zeroDepot, endDepot, N_type = build_node_id_sets(nodes)
# #     if use_imjn:
# #         R = build_requests(nodes, use_imjn=use_imjn)
# #         N, P, D, C, F, zeroDepot, endDepot= build_node_imjn_sets(nodes,max_visits_transfer=max_visits_transfer, ev_constraints=ev_constraints)
# #         N_minor, P_minor, D_minor, C_minor, F_minor, zeroDepot_bin, endDepot_bin= build_node_id_sets(nodes)
# #         Departures = build_departures(t_transfer, C_minor)
# #     else: 
# #         R, Cr = build_requests(nodes, use_imjn=use_imjn)
# #         N, P, D, C, F, zeroDepot, endDepot = build_node_id_sets(nodes, ev_constraints=ev_constraints)

# #     fi_r = build_fi_r(nodes, N, R, P, D)
# #     # print("fi_r is: ",fi_r)

# #     # Build parametres
# #     di, ei, li, qr, pair_pi_di, Lbar = build_service_params(nodes, t_transfer, P, D)

# #     # Constants
# #     Q = 4
# #     T = 3600 * 4
# #     V = n_vehicles
# #     K = list(range(V))
# #     ek = {k: 0 for k in K}
# #     lk = {k: 86400 for k in K}
# #     M = 100000
# #     n = len(P)
# #     # Electric Vehicle Constants
# #     C_bat_kWh = 100 # Battery Capacity in kWh
# #     C_bat = C_bat_kWh*(60*60) # Battery Capacity in kWs
# #     alpha = 200 # Charging Power in kW
# #     beta = 15 # Average Power in kW
# #     gamma = 0.1 # Minimum battery capacity set at 10%

# #     # --- Dataset with all variants included ---

# #     sets = dict(
# #         # === Node & Set Information ===
# #         nodes=nodes, N=N, P=P, D=D, C=C, F=F, R=R, K=K,
# #         zeroDepot=zeroDepot, endDepot=endDepot)
# #     params = dict(
# #         # === Parameters ===
# #         Q=Q, T=T, qr=qr, di=di, ei=ei, li=li,
# #         ek=ek, lk=lk, fi_r=fi_r, Lbar=Lbar, pair_pi_di=pair_pi_di, M=M, n=n)

# #     params['tij'], params['cij'] = build_tij(
# #         t, nodes, n_requests, n_vehicles, n_trans_nodes,
# #         include_charging=ev_constraints, duplicate=duplicate_transfers,
# #         eliminate_arcs=arc_elimination,
# #         di=di, ei=ei, li=li, Lbar=Lbar, P=P, D=D, N=N,
# #         zeroDepot=zeroDepot, endDepot=endDepot, Cr=Cr if duplicate_transfers else None,
# #         n=n, pair_pi_di=pair_pi_di)
    
# #     tij_keys = params['tij'].keys()
# #     sets["A"] = {(i, j) for i in N for j in N if (self.base(i), self.base(j)) in tij_keys}

# #     if duplicate_transfers:
# #         params["Cr"] = Cr 
# #     else: params["Cr"] = None
# #     if ev_constraints:
# #         params["C_bat"] = C_bat
# #         params["alpha"] = alpha
# #         params["beta"] = beta
# #         params["gamma"] = gamma
# #     if use_imjn:
# #         params['Departures'] = Departures
# #     else: params['Departures'] = None

# #     return sets, params


# # # ### === CHECKS === ###
# # sets1, params1 = pack_data()
# # print(params1['ei'])
# # print("-----")
# # print(params1['li'])
# # # print("fi_r is:", params1['fi_r'])
# # # print("This is the set A: ", sets1['A'])

# # # print("sets1 is :", sets1)

# # # print("Type of A:", type(sets1["A"]))
# # # if len(sets1["A"]) > 0:
# # #     first_el = next(iter(sets1["A"]))
# # #     print("First element of A:", first_el, "| Type:", type(first_el))

# # # sets2, params2 = pack_data(duplicate_transfers=False, use_imjn=True)
# # # print("N is :", sets2['N'])
# # # print("A is :", sets2['A'])
# # # print("Departures keys example:", list(params2["Departures"].keys())[:10])
# # # print("Types:", [type(k[0]) for k in params2["Departures"].keys()])
# # # print("nodes :", sets2["nodes"])
# # # print("C is: ", sets2["C"])
# # # print("Departures is: ", params2['Departures'].keys())

# # # sets3, params3 = pack_data(duplicate_transfers=False)

# # # def safe_get(d, key, default=None):
# # #     return d[key] if key in d else default

# # # print("\n========================")
# # # print("DATASET VALIDATION")
# # # print("========================\n")

# # # ### --- General Information --- ###
# # # print("=== General Node Info ===")
# # # print(f"Total nodes: {sets1['nodes'].shape[0]}")
# # # print(f"N={len(sets1['N'])}, P={len(sets1['P'])}, D={len(sets1['D'])}, "
# # #       f"C={len(sets1['C'])}, F={len(sets1['F'])}")
# # # print(f"ZeroDepot={sets1['zeroDepot']}, EndDepot={sets1['endDepot']}")
# # # print("-" * 40)

# # # ### --- Transfer / Duplication Check --- ###
# # # print("=== Transfer Duplication ===")
# # # if "Cr" in params1:
# # #     print("Cr dictionary present")
# # #     print("Cr keys:", list(params1["Cr"].keys()))
# # #     example_r = list(params1["Cr"].keys())[0]
# # #     print(f"Example Cr[{example_r}]:", params1["Cr"][example_r])
# # # else:
# # #     print("No Cr (duplicate_transfers=False)")
# # # print("-" * 40)

# # # ### --- Travel Time Dictionary --- ###
# # # print("=== Travel Time Dictionary ===")
# # # print(f"tij entries: {len(params1['tij'])}, cij entries: {len(params1['cij'])}")
# # # print("Example tij[(0,1)]:", params1["tij"].get((0, 1), "missing"))
# # # print(f"Max tij value: {max(params1['tij'].values()):.2f}")
# # # print(f"Min tij value: {min(params1['tij'].values()):.2f}")
# # # print("-" * 40)

# # # ### --- Request Consistency --- ###
# # # print("=== Request and Flow Checks ===")
# # # print("Requests R:", sets1["R"])
# # # print("Example pair_pi_di:", list(params1["pair_pi_di"].items())[:3])
# # # print("Example fi_r entries:", list(params1["fi_r"].items())[:5])
# # # print("-" * 40)

# # # ### --- EV Parameter Check --- ###
# # # print("=== EV Parameters ===")
# # # if "C_bat" in params1:
# # #     print(f"C_bat={params1['C_bat']}, α={params1['alpha']}, β={params1['beta']}, γ={params1['gamma']}")
# # # else:
# # #     print("No EV parameters (ev_constraints=False).")
# # # print("-" * 40)

# # # ### --- IMJN Departures (if applicable) --- ###
# # # print("=== IMJN Departures ===")
# # # if "Departures" in params2 and isinstance(params2["Departures"], dict):
# # #     dep = params2["Departures"]
# # #     print(f"Number of arcs in Departures: {len(dep)}")
# # #     first_arc, times = next(iter(dep.items()))
# # #     print("Example arc:", first_arc, "→", list(times.items())[:3])
# # # else:
# # #     print("No Departures (use_imjn=False).")
# # # print("-" * 40)

# # # ### --- Comparative Summary --- ###
# # # print("=== Comparative Summary ===")
# # # print(f"Default (dup=True): N={len(sets1['N'])}, C={len(sets1['C'])}")
# # # print(f"IMJN (use_imjn=True): N={len(sets2['N'])}, C={len(sets2['C'])}")
# # # print(f"No duplication: N={len(sets3['N'])}, C={len(sets3['C'])}")
# # # print("========================\n CHECK COMPLETE ")




# ########## Constraints ##########

# # # def base(i):
# # #     return i[0] if isinstance(i, tuple) else i


# # # # === Decision Variables ===
# # # def define_variables(m, sets, params, 
# # #     ev_constraints=False,
# # #     variable_substitution=True,
# # #     Timetabled_Departures=False,
# # #     use_imjn=True):

# # #     """
# # #     Define decision variables for the IDARP model, adapting to:
# # #       - imjn artificial nodes (i,m)
# # #       - timetabled departures (with z[d,i,j])
# # #       - EV constraints (battery and charging)
# # #       - variable substitution (link x and v)
# # #     """

# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]

# # #     Departures = params.get("Departures", None)

# # #     ### --- Original variables ---
# # #     # Vehicle routing arc variables
# # #     x = m.addVars(K, A, vtype=gb.GRB.BINARY, name="x")   # vehicle arcs
# # #     # Passenger flow
# # #     y = m.addVars(R, A, vtype=gb.GRB.CONTINUOUS, lb=0, name="y")
# # #     # Node timing variables
# # #     T_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="T_node")
# # #     T_veh  = m.addVars(K, vtype=gb.GRB.CONTINUOUS, name="T_veh")
# # #     D_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="D_node")

# # #     ### --- Developments ---
# # #     # Optional variables
# # #     v, a, B_node, E_node = None, None, None, None
# # #     z = {}

# # #     # Variable Substitution
# # #     if variable_substitution:
# # #         v = m.addVars(A, vtype=gb.GRB.BINARY, name="v")      # DAR arcs

# # #     # Transfer / fixed-route arcs
# # #     z = {}
# # #     if Timetabled_Departures and Departures is not None:
# # #         # Timetabled departures case
# # #         for (i, j) in A:
# # #             if i in C and j in C and (i[0], j[0]) in Departures:
# # #                 for d in Departures[(i[0], j[0])]:
# # #                     z[(d, i, j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{d},{i},{j}]")
# # #             # Artificial node marker (used in some strengthened formulations)

# # #     if use_imjn:
# # #         a = m.addVars(N, vtype=gb.GRB.BINARY, name="a")

# # #         # print("Checking problematic z keys...")
# # #         # if (1, (11,0), (10,0)) in z:
# # #         #     print("✅ Found z[(1,(11,0),(10,0))]")
# # #         # else:
# # #         #     print("❌ Missing z[(1,(11,0),(10,0))]")

# # #         # # Check other direction
# # #         # if (1, (10,0), (11,0)) in z:
# # #         #     print("✅ Found z[(1,(10,0),(11,0))]")
# # #         # else:
# # #         #     print("❌ Missing z[(1,(10,0),(11,0))]")

# # #         # # Check what d values exist
# # #         # for key in z.keys():
# # #         #     if key[1] == (11,0) and key[2] == (10,0):
# # #         #         print("Found matching z:", key)


# # #     elif C:  
# # #         # Transfer Original Modeling
# # #         for i in C:
# # #             for j in C:
# # #                 z[(i, j)] = m.addVar(vtype=gb.GRB.BINARY, name=f"z[{i},{j}]")

# # #     # EV-specific variables
# # #     if ev_constraints:
# # #         B_node = m.addVars(N, vtype=gb.GRB.CONTINUOUS, name="B_node")  # battery SOC
# # #         E_node = m.addVars(F, vtype=gb.GRB.CONTINUOUS, name="E_node")  # recharge time

# # #     return {
# # #         "x": x, "v": v, "y": y, "z": z,
# # #         "T_node": T_node, "T_veh": T_veh, "D_node": D_node,
# # #         "a": a, "B_node": B_node, "E_node": E_node
# # #     }


# # # # === Objective ===
# # # def set_objective(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     cij = params["cij"]
# # #     x, v = vars_['x'], vars_["v"]
# # #     # ride_time = 
# # #     if variable_substitution:
# # #         m.setObjective(
# # #             gb.quicksum(cij[(base(i),base(j))] * v[i,j] for (i,j) in A),
# # #             gb.GRB.MINIMIZE
# # #         )

# # #     else: m.setObjective(gb.quicksum(cij[i[0] if isinstance(i, tuple) else i, j[0] if isinstance(j, tuple) else j] * x[k,i,j] for k in K for (i,j) in A), gb.GRB.MINIMIZE)

# # # # === Constraint groups ===

# # # def add_vehicle_logic_constraints(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     v, x, y = vars_["v"], vars_["x"], vars_["y"]
# # #     zeroDepot_node, endDepot_node = zeroDepot, endDepot

# # #     if variable_substitution:
# # #         # pickups/dropoffs visited exactly once
# # #         for i in P + D:
# # #             m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in A) == 1)
# # #             m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) == 1)

# # #         # transfer nodes at most once
# # #         for i in C:
# # #             m.addConstr(gb.quicksum(v[i,j] for j in N if (i,j) in A) <= 1)
# # #             m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) <= 1)

# # #         # depot logic
# # #         for k in K:
# # #             m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for j in N if (zeroDepot_node,j) in A) <= 1)
# # #         for j in N:
# # #             if (zeroDepot_node, j) in A:
# # #                 m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for k in K) == v[zeroDepot_node,j])

# # #         # flow conservation
# # #         for i in (P+D+C):
# # #             m.addConstr(gb.quicksum(v[j,i] for j in N if (j,i) in A) -
# # #                         gb.quicksum(v[i,j] for j in N if (i,j) in A) == 0)

# # #         # vehicle must end at depot
# # #         m.addConstr(gb.quicksum(v[i,endDepot_node] for i in N if (i,endDepot_node) in A) == len(K))
# # #         for k in K:
# # #             m.addConstr(gb.quicksum(x[k,i,endDepot_node] for i in N if (i,endDepot_node) in A) == 1)

# # #         for i in N:
# # #             if (i,endDepot_node) in A:
# # #                 m.addConstr(
# # #                     gb.quicksum(x[k,i,endDepot_node] for k in K) == v[i,endDepot_node],
# # #                     name=f"endDepot_v_link[{i}]"
# # #                 )

# # #         # added by me
# # #         for r in R:
# # #             for (i, j) in A:
# # #                 if (r, i, j) in y:
# # #                     # y[r,i,j] ≤ v[i,j]
# # #                     m.addConstr(y[r, i, j] <= v[i, j],
# # #                                 name=f"passenger_vehicle_link[{r},{i},{j}]")

# # #     else: 
# # #         # pickups/dropoffs visited exactly once by some vehicle
# # #         for i in P + D:
# # #             m.addConstr(gb.quicksum(x[k,i,j] for k in K for j in N if (i,j) in A) == 1)

# # #         # each vehicle leaves depot exactly once
# # #         for k in K:
# # #             m.addConstr(gb.quicksum(x[k,zeroDepot_node,j] for j in N if (zeroDepot_node,j) in A) == 1)

# # #         # vehicle flow conservation
# # #         for k in K:
# # #             for i in (P + D + C):
# # #                 m.addConstr(gb.quicksum(x[k,j,i] for j in N if (j,i) in A) -
# # #                             gb.quicksum(x[k,i,j] for j in N if (i,j) in A) == 0)

# # #         # each vehicle ends at depot
# # #         for k in K:
# # #             m.addConstr(gb.quicksum(x[k,i,endDepot_node] for i in N if (i,endDepot_node) in A) == 1)

# # #         # added by me
# # #         for r in R:
# # #             for (i, j) in A:
# # #                 if (r, i, j) in y:
# # #                     # y[r,i,j] ≤ v[i,j]
# # #                     m.addConstr(y[r, i, j] <= gb.quicksum(x[k, i, j] for k in K),
# # #                                 name=f"passenger_vehicle_link[{r},{i},{j}]")



# # # def add_passenger_balance_constraints(m, vars_, sets, params, Timetabled_Departures=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     y, z = vars_["y"], vars_["z"]
# # #     Cr, fi_r, Departures = params["Cr"], params["fi_r"], params["Departures"]

# # #     if Timetabled_Departures and Departures != None:
# # #         for r in R:
# # #             for i in C:
# # #                 lhs_out = gb.quicksum(y[r, i, j] for j in N if (i, j) in A)
# # #                 lhs_in  = gb.quicksum(y[r, j, i] for j in N if (j, i) in A)
# # #                 z_out = gb.quicksum(z[(d, i, j)] for j in C if (i, j) in A and i != j and (base(i), base(j)) in Departures for d in Departures[(base(i), base(j))].keys())
# # #                 z_in  = gb.quicksum(z[(d, j, i)] for j in C if (i, j) in A and i != j and (base(j), base(i)) in Departures for d in Departures[(base(j), base(i))].keys())
# # #                 m.addConstr(lhs_out + z_out - lhs_in - z_in == fi_r[(r, i)], name=f"pass_bal_transfer[{r},{i}]")
# # #     else: 
# # #         if Cr:
# # #             for r in R:
# # #                 for i in Cr[r]:
# # #                     m.addConstr(
# # #                         gb.quicksum(y[r,i,j] for j in N if (i,j) in A) + gb.quicksum(z[(i,j)] for j in Cr[r] if (i,j) in A)
# # #                         - gb.quicksum(y[r,j,i] for j in N if (j,i) in A) - gb.quicksum(z[(j,i)] for j in Cr[r]if (j,i) in A)
# # #                         == fi_r[(r,i)],
# # #                         name=f"pass_bal_transfer[{r},{i}]"
# # #                     )

# # #     # non-transfer nodes
# # #     for r in R:
# # #         for i in N:
# # #             if i not in C:
# # #                 m.addConstr(
# # #                     gb.quicksum(y[r,i,j] for j in N if (i,j) in A) - gb.quicksum(y[r,j,i] for j in N if (j,i) in A) == fi_r[(r,i)],
# # #                     name=f"pass_bal_nontransfer[{r},{i}]"
# # #                 )


# # # def add_capacity_constraints(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     x, v, y = vars_['x'], vars_["v"], vars_["y"]
# # #     qr, Q = params["qr"], params["Q"]
# # #     if variable_substitution:
# # #         for (i,j) in A:
# # #             m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * v[i,j])

# # #     else:
# # #         for (i,j) in A:
# # #             m.addConstr(gb.quicksum(qr[r] * y[r,i,j] for r in R) <= Q * gb.quicksum(x[k,i,j] for k in K))


# # # def add_time_consistency_constraints(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     x, v, z, T_node, T_veh, D_node = vars_["x"], vars_["v"], vars_["z"], vars_["T_node"], vars_["T_veh"], vars_["D_node"]
# # #     tij, di, ei, li, ek, lk, Lbar, M, T, Cr = params["tij"], params["di"], params["ei"], params["li"], params["ek"], params["lk"], params['Lbar'], params["M"], params["T"], params['Cr']
# # #     pair_pi_di = params["pair_pi_di"]
# # #     zeroDepot_node, endDepot_node = zeroDepot, endDepot

# # #     if variable_substitution:
# # #         # depot to first node
# # #         for k in K:
# # #             for j in N:
# # #                 if (zeroDepot_node,j) in A:
# # #                     m.addConstr(T_node[j] >= T_veh[k] + tij[(base(zeroDepot_node), base(j))] - M*(1-x[k,zeroDepot_node,j]))

# # #         # service times
# # #         for (i,j) in A:
# # #             m.addConstr(T_node[j] >= T_node[i] + di[base(i)] + tij[(base(i),base(j))] - M*(1-v[i,j]))

# # #         # ride times
# # #         for p, d in pair_pi_di.items():
# # #             m.addConstr(T_node[d] - (T_node[p] + di[base(p)]) >= tij[(base(p), base(d))])
# # #             m.addConstr(T_node[d] - (T_node[p] + di[base(p)]) <= Lbar[base(p)])

# # #         # time windows
# # #         for i in P+D:
# # #             m.addConstr(T_node[i] >= ei[base(i)])
# # #             m.addConstr(T_node[i] <= li[base(i)])

# # #         # cumulative time
# # #         for (i,j) in A:
# # #             m.addConstr(D_node[j] >= D_node[i] + di[base(i)] + tij[(base(i),base(j))] - M*(1-v[i,j]))

# # #         # max route duration
# # #         for i in N:
# # #             if (i,endDepot_node) in A:
# # #                 m.addConstr(D_node[i] + di[base(i)] + tij[(base(i),base(endDepot_node))] <= T)

# # #         # depot departure window
# # #         for k in K:
# # #             m.addConstr(T_veh[k] >= ek[k])
# # #             m.addConstr(T_veh[k] <= lk[k])

# # #         for r in R:
# # #             if Cr:
# # #                 for i in Cr[r]:
# # #                     for j in Cr[r]:
# # #                         if (i,j) in tij:
# # #                             m.addConstr(
# # #                                 T_node[j] >= T_node[i] + tij[base(i),base(j)] - M*(1 - z[(i,j)]),
# # #                                 name=f"time_link_fixed[{r},{i},{j}]"
# # #                             )
   
# # #     else:
# # #         # (9) Timing from depot to first served node on vehicle k
# # #         for k in K:
# # #             for i in N:
# # #                 for j in N:
# # #                     if (i,j) in tij and (zeroDepot_node,j) in tij:
# # #                         m.addConstr(
# # #                             T_node[j] >= T_veh[k] + tij[base(zeroDepot_node),base(j)] - M*(1 - x[k,zeroDepot_node,j]),
# # #                             name=f"time_dep_first[{k},{j}]"
# # #                         )

# # #         # (10) Timing between consecutive service nodes (DRT)
# # #         for i in N:
# # #             for j in N:
# # #                 if (base(i),base(j)) in tij:
# # #                     m.addConstr(
# # #                         T_node[j] >= T_node[i] + di[base(i)] + tij[base(i),base(j)] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
# # #                         name=f"time_link_drt[{i},{j}]"
# # #                     )

# # #         # (11) Timing along fixed-route between transfer nodes for request r
# # #         for r in R:
# # #             if Cr:
# # #                 for i in Cr[r]:
# # #                     for j in Cr[r]:
# # #                         if (i,j) in A:
# # #                             m.addConstr(
# # #                                 T_node[j] >= T_node[i] + tij[base(i),base(j)] - M*(1 - z[(i,j)]),
# # #                                 name=f"time_link_fixed[{r},{i},{j}]"
# # #                             )

# # #         # (12) Time windows at pickup and dropoff nodes
# # #         ### for i in (P + D): ### No time windows for dropoffs in this example
# # #         for i in P + D:
# # #             m.addConstr(T_node[i] >= ei[base(i)], name=f"tw_lo[{i}]")
# # #             m.addConstr(T_node[i] <= li[base(i)], name=f"tw_hi[{i}]")

# # #         # (13) Ride-time constraints: t_{i,n+i} <= B_{n+i} - (B_i + d_i) <= Lbar_i  for pickups i
# # #         for i in P:
# # #             drop = pair_pi_di[i]  # dropoff node for pickup i (usually n+i)
# # #             if (i,drop) in A:
# # #                 m.addConstr(T_node[drop] - (T_node[i] + di[base(i)]) >= tij[base(i),base(drop)], name=f"ride_lo[{i}]")
# # #                 m.addConstr(T_node[drop] - (T_node[i] + di[base(i)]) <= Lbar[base(i)],   name=f"ride_hi[{i}]")

# # #         # (14) Cumulative time variable D (mirrors (10) without B, used for duration in (15))
# # #         for i in N:
# # #             for j in N:
# # #                 if (i,j) in A:
# # #                     m.addConstr(
# # #                         D_node[j] >= D_node[i] + di[base(i)] + tij[base(i),base(j)] - M*(1 - gb.quicksum(x[k,i,j] for k in K)),
# # #                         name=f"D_link[{i},{j}]"
# # #                     )

# # #         # (15) Vehicle maximum route duration T (ensure from node i to end depot fits in T)
# # #         for i in N:
# # #             if (i,endDepot_node) in A:
# # #                 m.addConstr(D_node[i] + di[base(i)] + tij[base(i),base(endDepot_node)] <= T, name=f"duration[{i}]")

# # #         # (16) Depot departure window for each vehicle
# # #         for k in K:
# # #             m.addConstr(T_veh[k] >= ek[k], name=f"veh_dep_lo[{k}]")
# # #             m.addConstr(T_veh[k] <= lk[k], name=f"veh_dep_hi[{k}]")


# # # def add_variable_substitution_constraints(m, vars_, sets, params):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     v, x = vars_["v"], vars_["x"]
# # #     zeroDepot_node = zeroDepot

# # #     for (i,j) in A:
# # #         if i != zeroDepot_node:
# # #             m.addConstr(v[i,j] == gb.quicksum(x[k,i,j] for k in K))


# # # def add_subtour_elimination_constraints(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     x, v = vars_['x'], vars_["v"]

# # #     if variable_substitution:
# # #         for i in N:
# # #             for j in N:
# # #                 if i != j and (i,j) in A and (j,i) in A:
# # #                     m.addConstr(v[i,j] + v[j,i] <= 1)
# # #     else:
# # #         for i in N:
# # #             for j in N:
# # #                 if i != j and (i,j) in A and (j,i) in A:
# # #                     m.addConstr(gb.quicksum(x[k,i,j] + x[k,j,i] for k in K) <= 1)


# # # def add_transfer_node_constraints(m, vars_, sets, params, variable_substitution=True, Timetabled_Departures=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     x, v, z = vars_['x'], vars_["v"], vars_["z"]
# # #     Cr, Departures = params['Cr'], params["Departures"]

# # #     if variable_substitution:
# # #         if Timetabled_Departures: 
# # #             for i in C:
# # #                 m.addConstr(
# # #                     gb.quicksum(v[j,i] for j in N if (j,i) in A) ==
# # #                     gb.quicksum(z[(d,i,j)] + z[(d,j,i)] for j in C if (i,j) in A and i!=j and (base(i),base(j)) in Departures for d in Departures[(base(i),base(j))].keys())
# # #                 )
# # #         else:
# # #             if Cr:
# # #                 for r in R:
# # #                     for i in Cr[r]:
# # #                         m.addConstr(
# # #                             gb.quicksum(v[j,i] for j in N if (j,i) in A) ==
# # #                             gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) in A)
# # #                         )
# # #                 for r in R:
# # #                     m.addConstr(gb.quicksum(v[i,j] for i in Cr[r] for j in N if (i,j) in A) <= 2)
# # #     else:
# # #         if Cr:
# # #             for r in R:
# # #                 for i in Cr[r]:
# # #                     m.addConstr(
# # #                         gb.quicksum(gb.quicksum(x[k,j,i] for k in K) for j in N if (j,i) in A) ==
# # #                         gb.quicksum(z[(i,j)] + z[(j,i)] for j in Cr[r] if (i,j) in A)
# # #                     )
# # #         if Cr:
# # #             for r in R:
# # #                 m.addConstr(gb.quicksum(gb.quicksum(x[k,i,j] for k in K) for i in Cr[r] for j in N if (i,j) in A) <= 2)    
    


# # # def add_battery_constraints(m, vars_, sets, params, variable_substitution=True):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     x, v, B_node, E_node = vars_['x'], vars_["v"], vars_["B_node"], vars_["E_node"]

# # #     tij = params["tij"]
# # #     beta, alpha, C_bat, gamma = params["beta"], params["alpha"], params["C_bat"], params["gamma"]

# # #     if variable_substitution:
# # #         for (i, j) in A:
# # #             if i not in F and j not in F:
# # #                 m.addConstr(B_node[j] - B_node[i] + beta * tij[(base(i), base(j))] - C_bat * (1 - v[i, j]) <= 0)
# # #                 m.addConstr(B_node[j] - B_node[i] + beta * tij[(base(i), base(j))] + C_bat * (1 - v[i, j]) >= 0)
# # #     else:
# # #         for (i, j) in A:
# # #             if i not in F and j not in F:
# # #                 m.addConstr(B_node[j] - B_node[i] + beta * tij[(base(i), base(j))] - C_bat * (1 - gb.quicksum(x[k,i,j] for k in K)) <= 0)

# # #     # battery capacity
# # #     for s in F:
# # #         m.addConstr(C_bat - (B_node[s] + alpha * E_node[s]) >= 0)

# # #     # min SOC
# # #     for i in N:
# # #         m.addConstr(B_node[i] >= gamma * C_bat)

# # #     # init SOC
# # #     m.addConstr(B_node[zeroDepot] == C_bat)


# # # def add_scheduled_PT_constraints(m, vars_, sets, params):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     z, T_node = vars_["z"], vars_["T_node"]
# # #     M, Departures, tij = params["M"], params["Departures"], params['tij']

# # #     for i in N:
# # #         if i not in C:
# # #             m.addConstr(
# # #                 M*(1 - gb.quicksum(z[(d,i,j)] for j in C if (i, j) in A and i != j and (base(j), base(i)) in Departures for d in Departures[(base(i),base(j))].keys())) +
# # #                 gb.quicksum(Departures[(base(i),base(j))][d] * z[(d,i,j)] for j in C if (i,j) in A and i!=j and (base(i),base(j)) in Departures for d in Departures[(base(i),base(j))].keys())
# # #                 - T_node[i] >= 0
# # #             )
# # #     for i in C:
# # #         expr = gb.quicksum(
# # #             (Departures[(base(i), base(j))][d] + tij[base(i), base(j)]) * z[d, i, j]
# # #             for j in C
# # #             if (i, j) in A
# # #             and (base(i), base(j)) in Departures
# # #             for d in Departures[(base(i), base(j))].keys()
# # #         )
# # #         m.addConstr(T_node[i] - expr >= 0, name=f"PT_departure_after_arrival_{i}")



# # # def add_artificial_node_constraints(m, vars_, sets, params):
# # #     nodes, N, P, D, C, F, R, K, zeroDepot, endDepot, A = sets["nodes"], sets["N"], sets["P"], sets["D"], sets["C"], sets["F"], sets["R"], sets["K"], sets["zeroDepot"], sets["endDepot"], sets["A"]
# # #     a = vars_["a"]
# # #     for (i,m_) in N:
# # #         if m_ < max(n for (j,n) in N if j == i):
# # #             m.addConstr(a[i,m_] - a[i,m_+1] >= 0)

# ######### Model ##########


# # # def build_model(
# # #     duplicate_transfers=False,
# # #     arc_elimination=False, 
# # #     variable_substitution=True,
# # #     subtour_elimination=True,
# # #     transfer_node_strengthening=False,
# # #     ev_constraints=False,
# # #     Timetabled_Departures=False,
# # #     use_imjn=False,
# # #     model_name="IDARP_Model",
# # #     TIME_LIMIT=10, #2*60*60,
# # #     sets={}, params={}
# # # ):
# # #     """
# # #     Build a Gurobi model under different structural conditions.
# # #     """

# # #     # # === Step 1: Unpack Data into (sets, params) ===
# # #     # sets, params = pack_data(
# # #     #     duplicate_transfers=duplicate_transfers,
# # #     #     arc_elimination=arc_elimination,
# # #     #     ev_constraints=ev_constraints,
# # #     #     use_imjn=use_imjn
# # #     # )

# # #     # === Step 2: Build Model ===
# # #     m = gb.Model(model_name)
# # #     m.setParam("TimeLimit", TIME_LIMIT)

# # #     # === Step 3: Define Variables ===
# # #     vars_ = define_variables(
# # #         m, sets, params,
# # #         ev_constraints=ev_constraints,
# # #         variable_substitution=variable_substitution,
# # #         Timetabled_Departures=Timetabled_Departures,
# # #         use_imjn=use_imjn
# # #     )
# # #     # # ======== Known routes ========
# # #     # known_vehicle_routes = {
# # #     #     1: [0, 2, 6, 3, 7, 9],
# # #     #     2: [0, 1, 5, 4, 8, 9],
# # #     # }

# # #     # known_request_routes = {
# # #     #     1: [1, 5],
# # #     #     2: [2, 6],
# # #     #     3: [3, 7],
# # #     #     4: [4, 8],
# # #     # }

# # #     # # ======== Initialize Starts ========
# # #     # v = vars_["v"]
# # #     # y = vars_["y"]

# # #     # # Reset all starts to 0
# # #     # for (i, j) in v.keys():
# # #     #     v[(i, j)].Start = 0
# # #     # for (r, i, j) in y.keys():
# # #     #     y[(r, i, j)].Start = 0

# # #     # # ======== Activate vehicle arcs ========
# # #     # for route in known_vehicle_routes.values():
# # #     #     for a, b in zip(route[:-1], route[1:]):
# # #     #         if (a, b) in v:
# # #     #             v[(a, b)].Start = 1

# # #     # # ======== Activate request arcs ========
# # #     # for r, route in known_request_routes.items():
# # #     #     for a, b in zip(route[:-1], route[1:]):
# # #     #         if (r, a, b) in y:
# # #     #             y[(r, a, b)].Start = 1

# # #     # print("✅ Warm start: loaded vehicle and request routes.")


# # #     # === Step 4: Objective ===
# # #     set_objective(m, vars_, sets, params, self.variable_substitution=variable_substitution)

# # #     # === Step 5: Constraints ===
# # #     add_vehicle_logic_constraints(m, vars_, sets, params, variable_substitution=variable_substitution)
# # #     add_passenger_balance_constraints(m, vars_, sets, params, Timetabled_Departures=Timetabled_Departures)
# # #     add_capacity_constraints(m, vars_, sets, params, variable_substitution=variable_substitution)
# # #     add_time_consistency_constraints(m, vars_, sets, params, variable_substitution=variable_substitution)

# # #     if variable_substitution:
# # #         add_variable_substitution_constraints(m, vars_, sets, params)

# # #     if subtour_elimination:
# # #         add_subtour_elimination_constraints(m, vars_, sets, params, variable_substitution)

# # #     if transfer_node_strengthening:
# # #         if duplicate_transfers or Timetabled_Departures:
# # #             add_transfer_node_constraints(m, vars_, sets, params,
# # #                                         variable_substitution=variable_substitution,
# # #                                         Timetabled_Departures=Timetabled_Departures)

# # #     if ev_constraints:
# # #         add_battery_constraints(m, vars_, sets, params, variable_substitution=variable_substitution)

# # #     if Timetabled_Departures:
# # #         add_scheduled_PT_constraints(m, vars_, sets, params)

# # #     if use_imjn:
# # #         add_artificial_node_constraints(m, vars_, sets, params)

# # #     # === Finalize ===
# # #     m.update()
# # #     return m, vars_  #, sets, params


# ############# Runner_full ##############
# ##### Imports #####
# from tqdm import tqdm
# import time as t
# import itertools
# import gurobipy as gb
# import pandas as pd

# from Model import build_model
# from routes_2 import extract_route_final
# from Parametres import pack_data
# from Debugger_code import save_model_debug_data, save_constraint_lhs_rhs
# from Cluster_Heuristic import subtour_callback

# TIME_LIMIT = 2 * 60 * 60  # 2 hours
# # ===================== CHOOSE PARAMETERS =====================
# # Set each parameter manually (True/False)
# bool_params_singular = {
#     "duplicate_transfers": False,
#     "arc_elimination": True,
#     "variable_substitution": True,
#     "subtour_elimination": True,
#     "transfer_node_strengthening": True,
#     "ev_constraints": False,
#     "timetabled_departures": True,
#     "use_imjn": True,
# }

# def run_singular_model(TIME_LIMIT, bool_params_singular):
#     # Model name for identification
#     model_name = "IDARP_Model_Single"

#     print(f"\n=== Running single model: {model_name} ===")
#     print(bool_params_singular)

#     # Impossible combinations


#     # Get data for case
#     data_sets, data_params = pack_data(duplicate_transfers=bool_params_singular["duplicate_transfers"], 
#                              arc_elimination=bool_params_singular["arc_elimination"], 
#                              ev_constraints=bool_params_singular["ev_constraints"], 
#                              use_imjn=bool_params_singular["use_imjn"])

#     # ===================== BUILD AND SOLVE =====================
#     start_cpu = t.process_time()
#     start_wall = t.perf_counter()

#     # Build & solve model
#     m, vars_ = build_model(
#         model_name=model_name,
#         TIME_LIMIT=TIME_LIMIT,
#         sets=data_sets,
#         params=data_params,
#         **bool_params_singular
#         )
    
#     # Attach required data for the callback
#     m._v = vars_["v"]              # aggregated arc variables
#     m._x = vars_["x"]              # per-vehicle arc variables
#     m._P = data_sets["P"]          # pickup nodes
#     m._D = data_sets["D"]          # pickup nodes
#     m._N = data_sets["N"]          # all nodes
#     m._K = data_sets["K"]          # vehicle set
#     m._A = data_sets["A"]          # feasible arcs
#     m.Params.LazyConstraints = 1   # enable lazy cuts


#     if bool_params_singular['subtour_elimination']:
#         m.optimize(lambda model, where: subtour_callback(model, where, variable_substitution=bool_params_singular['variable_substitution']))
#     else:
#         m.optimize()

#     # if m.SolCount > 0:
#     #     print(f"✅ Feasible solution found (Obj = {m.ObjVal:.2f})")
#     #     save_model_debug_data(m, vars_, model_name=model_name)
#     #     save_constraint_lhs_rhs(m, file_name=f"{model_name}_constraints.csv")
#     # else:
#     #     print("⚠️ No feasible solution found or time limit reached — skipping LHS/RHS save.")

#     end_cpu = t.process_time()
#     end_wall = t.perf_counter()
#     cpu_time = end_cpu - start_cpu
#     elapsed_time = end_wall - start_wall

#     # ===================== HANDLE INFEASIBILITY =====================
#     if m.status == gb.GRB.INFEASIBLE:
#         print("⚠️ Model infeasible — writing IIS files...")
#         m.computeIIS()
#         m.write(f"{model_name}.ilp")
#         m.write(f"{model_name}.iis")
#         raise SystemExit("Model infeasible; IIS written.")

#     # ===================== EXTRACT ROUTES =====================
#     v1, v2, r1, r2, r3, r4 = extract_route_final(
#         vars_,
#         nodes=data_sets["nodes"],
#         N=data_sets["N"],
#         K=data_sets["K"],
#         variable_substitution=bool_params_singular['variable_substitution']
#     )

#     # ===================== SAVE RESULTS =====================
#     results = [{
#         "Model": model_name,
#         **bool_params_singular,
#         "Objective": m.ObjVal,
#         "Lower Bound": m.ObjBound,
#         "Vehicle 1 Route": str(v1),
#         "Vehicle 2 Route": str(v2),
#         "Request 1 Route": str(r1),
#         "Request 2 Route": str(r2),
#         "Request 3 Route": str(r3),
#         "Request 4 Route": str(r4),
#         "CPU Time (s)": cpu_time,
#         "Elapsed Time (s)": elapsed_time,
#     }]
    
#     print(results)

#     # df = pd.DataFrame(results)
#     # df.to_csv("results_single.csv", index=False)

#     # print(f"\n✅ Done: {model_name} | Obj = {m.ObjVal:.2f}, Bound = {m.ObjBound:.2f}")
#     # print("✅ Results saved to results_single.csv")


# def run_all_possibilities(TIME_LIMIT):
#     # ===================== PARAMETERS =====================
#     boolean_params__multiple = {
#         'duplicate_transfers': [True, False],
#         'arc_elimination': [True, False],
#         'variable_substitution': [True, False],
#         'subtour_elimination': [True, False],
#         'timetabled_departures': [True, False],
#         'use_imjn': [True, False],
#     }

#     results = []
#     skipped_combos = []

#     # Build all parameter combinations first
#     all_combos = list(itertools.product(*boolean_params__multiple.values()))
#     total_combos = len(all_combos)

#     print(f"=== Starting {total_combos} model configurations ===")

#     # Progress bar setup
#     pbar = tqdm(total=total_combos, desc="Progress", unit="config", ncols=100)
#     start_time_global = t.perf_counter()

#     # ===================== RUN ALL COMBINATIONS =====================
#     for idx, combo in enumerate(all_combos, start=1):
#         bool_params = dict(zip(boolean_params__multiple.keys(), combo))
#         model_name = "IDARP_Model_" + "_".join([f"{k[:3]}{int(v)}" for k, v in bool_params.items()])

#         # --- Default modifiers ---
#         bool_params['ev_constraints'] = False
#         bool_params['transfer_node_strengthening'] = True

#         # --- Feasibility filter ---
#         skip_reason = None
#         if bool_params["timetabled_departures"] and not bool_params["use_imjn"]:
#             skip_reason = "timetabled departures require IMJN nodes (use_imjn=True)."
#         elif bool_params["use_imjn"] and bool_params["duplicate_transfers"]:
#             skip_reason = "Duplicate transfers cannot be used with IMJN nodes."
#         elif bool_params["transfer_node_strengthening"] and not (
#             bool_params["timetabled_departures"] or bool_params["duplicate_transfers"]
#         ):
#             skip_reason = "Transfer node strengthening requires either Timetabled Departures or duplicate transfers."
#         elif bool_params["timetabled_departures"] and not bool_params["variable_substitution"]:
#             skip_reason = "Timetabled Departures without variable substitution are inconsistent."
#         elif not bool_params["duplicate_transfers"] and not bool_params["use_imjn"]:
#             skip_reason = "Cannot disable duplicate_transfers without IMJN nodes."

#         # --- Skip invalid combos ---
#         if skip_reason:
#             skipped_combos.append({**bool_params, "Reason": skip_reason})
#             pbar.update(1)
#             pbar.set_postfix({"Status": "⏭️ skipped"})
#             continue

#         # --- Start timer for ETA ---
#         start_wall = t.perf_counter()

#         # ===================== BUILD MODEL =====================
#         data_sets, data_params = pack_data(
#             duplicate_transfers=bool_params["duplicate_transfers"],
#             arc_elimination=bool_params["arc_elimination"],
#             ev_constraints=bool_params['ev_constraints'],
#             use_imjn=bool_params["use_imjn"]
#         )

#         m, vars_ = build_model(
#             model_name=model_name,
#             TIME_LIMIT=TIME_LIMIT,
#             sets=data_sets,
#             params=data_params,
#             **bool_params,
#         )

#         heuristic = DARPHeuristic(m, vars_, self.sets, self.params,
#                           variable_substitution=variable_substitution)

#         # Attach required data for the callback
#         m._v = vars_["v"]
#         m._x = vars_["x"]
#         m._P = data_sets["P"]
#         m._D = data_sets["D"]
#         m._N = data_sets["N"]
#         m._K = data_sets["K"]
#         m._A = data_sets["A"]
#         m.Params.LazyConstraints = 1

#         # ===================== SOLVE =====================
#         if bool_params['subtour_elimination']:
#             m.optimize(lambda model, where: heuristic.subtour_callback(model, where,
#                 variable_substitution=bool_params['variable_substitution']))
#         else:
#             m.optimize()

#         # ===================== RECORD RESULT =====================
#         end_wall = t.perf_counter()
#         elapsed_time = end_wall - start_wall

#         if m.status == gb.GRB.INFEASIBLE:
#             print(f"❌ {model_name} infeasible — skipping.")
#             results.append({
#                 "Model": model_name,
#                 **bool_params,
#                 "Objective": "Infeasible",
#                 "Lower Bound": "N/A",
#                 "Elapsed Time (s)": elapsed_time,
#             })
#         elif m.SolCount == 0:
#             print(f"⚠️ {model_name} — No feasible solution or timeout.")
#             results.append({
#                 "Model": model_name,
#                 **bool_params,
#                 "Objective": "No feasible solution",
#                 "Lower Bound": "N/A",
#                 "Elapsed Time (s)": elapsed_time,
#             })
#         else:
#             print(f"✅ {model_name} solved — Obj = {m.ObjVal:.2f}, Bound = {m.ObjBound:.2f}")
#             try:
#                 v1, v2, r1, r2, r3, r4 = extract_route_final(
#                     vars_,
#                     nodes=data_sets['nodes'],
#                     N=data_sets['N'],
#                     K=data_sets['K'],
#                     variable_substitution=bool_params["variable_substitution"]
#                 )
#             except Exception as e:
#                 print(f"⚠️ Route extraction failed: {e}")
#                 v1 = v2 = r1 = r2 = r3 = r4 = []

#             results.append({
#                 "Model": model_name,
#                 **bool_params,
#                 "Objective": m.ObjVal,
#                 "Lower Bound": m.ObjBound,
#                 "Vehicle 1 Route": str(v1),
#                 "Vehicle 2 Route": str(v2),
#                 "Request 1 Route": str(r1),
#                 "Request 2 Route": str(r2),
#                 "Request 3 Route": str(r3),
#                 "Request 4 Route": str(r4),
#                 "Elapsed Time (s)": elapsed_time,
#             })

#         # --- Update progress bar + ETA ---
#         avg_time = (t.perf_counter() - start_time_global) / idx
#         remaining = avg_time * (total_combos - idx)
#         eta_min, eta_sec = divmod(int(remaining), 60)
#         pbar.update(1)
#         pbar.set_postfix({
#             "ETA": f"{eta_min:02d}:{eta_sec:02d}",
#             "Last": f"{elapsed_time:.1f}s"
#         })

#     pbar.close()

#     # ===================== SAVE RESULTS =====================
#     df = pd.DataFrame(results)
#     df.to_csv("results_full.csv", index=False)
#     print("\n✅ All feasible results saved to results_full.csv")

#     if skipped_combos:
#         df_skipped = pd.DataFrame(skipped_combos)
#         df_skipped.to_csv("skipped_combinations.csv", index=False)
#         print(f"⚙️ Skipped {len(skipped_combos)} invalid configurations → saved to skipped_combinations.csv")


# run_singular_model(TIME_LIMIT, bool_params_singular)
# # run_all_possibilities(TIME_LIMIT)


