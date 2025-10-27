##### Imports #####
import numpy as np

class DARPDataBuilder:
    def __init__(self, duplicate_transfers=True, arc_elimination=True,
                 ev_constraints=False, use_imjn=False, MoPS = False, max_visits_transfer=5):
        self.duplicate_transfers = duplicate_transfers
        self.arc_elimination = arc_elimination
        self.ev_constraints = ev_constraints
        self.use_imjn = use_imjn
        self.MoPS = MoPS
        self.max_visits_transfer = max_visits_transfer

    def base(self, x):
        return x[0] if isinstance(x, tuple) else x
    
    def build_t_transfer(self, t, n_requests, n_vehicles, n_trans_nodes, nodes):
        """
        Build the extended travel-time matrix including transfer and/or charging nodes.

        Parameters
        ----------
        t : np.ndarray
            Base travel time matrix (n_base × n_base).
        n_requests : int
            Number of requests.
        n_vehicles : int
            Number of vehicles.
        n_trans_nodes : int
            Number of transfer nodes per cluster.
        duplicate_transfers : bool, optional
            If True, each request gets its own duplicated set of transfer nodes.
            If False, transfers are shared across requests.
        ev_constraints : bool, optional
            If True, also add charging station nodes for each vehicle.

        Returns
        -------
        t_transfer : np.ndarray
            The full expanded travel-time matrix including transfers and (optionally) charging stations.
        """

        n_base = t.shape[0]                # Base node count (P + D + depots, etc.)
        # print("n_base is: ",n_base)
        n_transfers = n_trans_nodes * ((n_requests-1) if self.duplicate_transfers else 0)
        # print("n_transfers:", n_transfers)
        n_charging = n_vehicles * n_trans_nodes if self.ev_constraints else 0
        # print("n_charging :", n_charging)
        n_total = n_base + n_transfers + n_charging
        # print("n_total", n_total)

        t_transfer = np.zeros((n_total, n_total))

        # ---- Base block ----
        t_transfer[:n_base, :n_base] = t
        # print("t_transfer[:n_base, :n_base]", t_transfer[:n_base, :n_base])

        # Index ranges for clarity
        # trans_start = n_base
        # trans_end = trans_start + n_transfers
        # charge_start = trans_end
        # charge_end = charge_start + n_charging

        dimension_base_matrix = max(set(nodes[:, 0])) + 1
        # print("dimension_base_matrix: ", dimension_base_matrix)

        # ---- Base → Transfers ----
        if n_transfers > 0:
            if self.ev_constraints and self.duplicate_transfers:
                n_columns = n_requests + n_vehicles - 1 ### Duplicate transfer nodes for each request and charging facilities for each vehicle
            elif self.ev_constraints: 
                n_columns = 1 ### Build transfer nodes and and duplicate them to have a charging facility at each transfer node
            elif self.duplicate_transfers:
                n_columns = n_requests - 1 ### duplicate alld transfer nodes for each request (no charging here)
            else: 
                n_columns = 0 ### Build transfer_nodes
                
            # print("n_columns", n_columns)

            base_to_trans = np.tile(t[:, dimension_base_matrix:dimension_base_matrix + n_columns],
                            (1, n_columns))
            # print("base_to_trans : ", base_to_trans)
            # print("t_transfer[:n_base, trans_start:trans_end]:", t_transfer[:n_base, n_base:n_base + n_transfers * n_columns])
            t_transfer[:n_base, n_base : n_base + n_transfers * n_columns] = base_to_trans

            trans_to_base = np.tile(t[dimension_base_matrix:dimension_base_matrix + n_columns, :],
                                    (n_columns, 1))
            t_transfer[n_base : n_base + n_transfers * n_columns, :n_base] = trans_to_base

            trans_block = np.tile(t[dimension_base_matrix:dimension_base_matrix + n_columns, dimension_base_matrix:dimension_base_matrix + n_columns],
                                (n_columns, n_columns))
            t_transfer[n_base : n_base + n_transfers * n_columns, n_base : n_base + n_transfers * n_columns] = trans_block

        # # ---- Base ↔ Charging stations ----
        # if self.ev_constraints and n_charging > 0:
        #     base_to_charge = np.tile(t[:, dimension_base_matrix:dimension_base_matrix + n_trans_nodes], (1, n_vehicles))
        #     charge_to_base = np.tile(t[dimension_base_matrix:dimension_base_matrix + n_trans_nodes, :], (n_vehicles, 1))
        #     charge_block = np.tile(t[dimension_base_matrix:dimension_base_matrix + n_trans_nodes, dimension_base_matrix:dimension_base_matrix + n_trans_nodes],
        #                         (n_vehicles, n_vehicles))

        #     # Base → Chargers
        #     t_transfer[:n_base, charge_start:charge_end] = base_to_charge
        #     # Chargers → Base
        #     t_transfer[charge_start:charge_end, :n_base] = charge_to_base
        #     # Chargers ↔ Chargers
        #     t_transfer[charge_start:charge_end, charge_start:charge_end] = charge_block

            # Optionally connect transfers ↔ chargers (only if both exist)
            # if n_transfers > 0:
            #     t_transfer[trans_start:trans_end, charge_start:charge_end] = np.mean(charge_block)
            #     t_transfer[charge_start:charge_end, trans_start:trans_end] = np.mean(charge_block)
        return t_transfer
    
    def build_nodes(self, base_nodes, n_requests, n_vehicles, n_trans_nodes):

        """
        Build node list with consistent numbering:
        - Transfer nodes start at 10, 11, 12, ...
        - Charging station nodes start right after the last transfer node.

        Parameters
        ----------
        base_nodes : np.ndarray
            Array of base nodes (depot, pickups, dropoffs)
        n_requests : int
            Number of requests
        n_vehicles : int
            Number of vehicles
        n_trans_nodes : int
            Number of physical transfer nodes per request or cluster
        duplicate_nodes : bool
            Whether to duplicate transfers for each request
        ev_constraints : bool
            Whether to include charging station nodes

        Returns
        -------
        nodes : np.ndarray
            Full node array including base, transfer, and optional charging nodes
        """
        transfers, charging_stations = [], []

        # --- TRANSFER NODES ---
        transfer_start_id = max(set(base_nodes[:, 0])) + 1
        if self.duplicate_transfers:
            # One set of transfers per request
            for r in range(1, n_requests + 1):
                for t_id in range(n_trans_nodes):
                    node_id = transfer_start_id + (r - 1) * n_trans_nodes + t_id
                    transfers.append([node_id, "transfer", str(r), "", "", 5, 0])
        else:
            # Single shared set of transfers
            for t_id in range(n_trans_nodes):
                node_id = transfer_start_id + t_id
                transfers.append([node_id, "transfer", "", "", "", 5, 0])

        # --- CHARGING STATIONS ---
        if self.ev_constraints:
            charge_start_id = transfer_start_id + n_trans_nodes * (n_requests if self.duplicate_transfers else 1)
            if self.duplicate_transfers:
                for k in range(1, n_vehicles + 1):
                    for c_id in range(n_trans_nodes):
                        node_id = charge_start_id + (k - 1) * n_trans_nodes + c_id
                        charging_stations.append([node_id, "charging station", "", 0, 86400, 0, 0])
            else:
                for c_id in range(n_trans_nodes):
                    node_id = charge_start_id + c_id
                    charging_stations.append([node_id, "charging station", "", 0, 86400, 0, 0])

        # --- SAFE CONVERSION ---
        if len(transfers) == 0:
            transfers = np.empty((0, base_nodes.shape[1]), dtype=object)
        if len(charging_stations) == 0:
            charging_stations = np.empty((0, base_nodes.shape[1]), dtype=object)

        # --- FINAL STACK ---
        return np.vstack([
            base_nodes,
            np.array(transfers, dtype=object),
            np.array(charging_stations, dtype=object)
        ])
  
    def build_node_id_sets(self, nodes):
        N = nodes[:, 0].astype(int).tolist()
        P = nodes[nodes[:, 1] == "pickup", 0].astype(int).tolist()
        D = nodes[nodes[:, 1] == "dropoff", 0].astype(int).tolist()
        C = nodes[nodes[:, 1] == "transfer", 0].astype(int).tolist()
        zeroDepot = int(nodes[nodes[:, 1] == "depot0", 0][0])
        endDepot = int(nodes[nodes[:, 1] == "depot1", 0][0])
        F = []
        P_M = []
        D_M = []
        if self.ev_constraints:
            F = nodes[nodes[:, 1] == "charging station", 0].astype(int).tolist()
        if self.MoPS:
            P_M = nodes[nodes[:, 1] == "MoPS pickup", 0].astype(int).tolist()
            D_M = nodes[nodes[:, 1] == "MoPS dropoff", 0].astype(int).tolist()
        return N, P, D, C, F, zeroDepot, endDepot, P_M, D_M

    def build_node_imjn_sets(self, nodes):
        N_minor, P_minor, D_minor, C_minor, F_minor, zeroDepot_node, endDepot_node, P_M_minor, D_M_minor = self.build_node_id_sets(nodes)
        # Wrap depots as (id,1)
        zeroDepot_node = (zeroDepot_node, 1)
        endDepot_node = (endDepot_node, 1)

        # Expand pickups/dropoffs as (i,1)
        P = [(i, 1) for i in P_minor]
        D = [(i, 1) for i in D_minor]

        # Expand transfer nodes with artificial visits
        C = [(i, m) for m in range(1, self.max_visits_transfer) for i in C_minor]
        F = []
        P_M = []
        D_M = []
        if self.ev_constraints:
            F = [(i, m) for i in (F_minor or []) for m in range(1, self.max_visits_transfer)]

        if self.MoPS:
            P_M = [(i, 1) for i in P_M_minor]
            D_M = [(i, 1) for i in D_M_minor]

        # All nodes in tuple format
        N = [zeroDepot_node, endDepot_node] + P + D + C + F + P_M + D_M

        return N, P, D, C, F, zeroDepot_node, endDepot_node, P_M, D_M

    def build_requests(self, nodes):
        if self.use_imjn:
            r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "transfer", "MoPS pickup", "MoPS dropoff"]), 2]
            r_values = r_values[r_values != ""]
            R = sorted(set(r_values.astype(int)))
            return R
        else:
            r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "transfer", "MoPS pickup", "MoPS dropoff"]), 2]
            r_values = r_values[r_values != ""]
            R = sorted(set(r_values.astype(int)))
            Cr = {r: nodes[(nodes[:, 1] == "transfer") & (nodes[:, 2] == str(r)), 0].astype(int).tolist()
                for r in R}
            return R, Cr

    def build_service_params(self, nodes, t, P, D, P_M, D_M):
        di = {int(row[0]): float(row[5]) for row in nodes if row[5] != ""}
        ei = {int(row[0]): float(row[3]) for row in nodes if row[3] != ""}
        li = {int(row[0]): float(row[4]) for row in nodes if row[4] != ""}
        r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "MoPS pickup", "MoPS dropoff"]), 2]
        r_values = r_values[r_values != ""]
        R = sorted(set(r_values.astype(int)))
        qr = {int(r): 1 for r in R}
        tij = {(i, j): float(t[i, j]) for i in range(len(t)) for j in range(len(t))}
        Lbar = {self.base(i): 1800 for i in P + P_M}
        pair_pi_di = {p: d for p, d in zip(P, D)}
        for p in pair_pi_di.keys():
            d = pair_pi_di[p]
            pid = self.base(p)
            did = self.base(d)
            e = ei[pid] + di[pid] + tij[pid, did]
            l = li[pid] + Lbar[pid]
            ei[did] = e
            li[did] = l
        pair_pi_di_M = {}
        if self.MoPS:
            pair_pi_di_M = {p: d for p, d in zip(P_M, D_M)}
            for p in pair_pi_di_M.keys():
                d = pair_pi_di_M[p]
                pid = self.base(p)
                did = self.base(d)
                e = ei[pid] + di[pid] + tij[pid, did]
                l = li[pid] + Lbar[pid]
                ei[did] = e
                li[did] = l
        return di, ei, li, qr, pair_pi_di, Lbar, pair_pi_di_M

    def build_fi_r(self, nodes, N, R, P, D, P_M, D_M):
        fi_r = {(int(r), i): 0 for r in R for i in N}
        for p, d in zip(P, D):
            r = int(nodes[nodes[:, 0] == self.base(p), 2][0])
            fi_r[r, p] = 1
            fi_r[r, d] = -1
        if self.MoPS:
            for p, d in zip(P_M, D_M):
                r = int(nodes[nodes[:, 0] == self.base(p), 2][0])
                fi_r[r, p] = 1
                fi_r[r, d] = -1
        return fi_r

    def build_departures(self, t, C_minor, interval=20, planning_horizon=24*60):    ##### Modify again ######
        """
        Build Departures dictionary with forward and backward arcs (including skip arcs).

        Args:
            t (np.array): travel time matrix (base or expanded).
            C_minor (list): list of transfer node ids (integers).
            interval (int): departure interval in seconds (default 300 = 5 minutes).
            planning_horizon (int): total time horizon (default 24h).

        Returns:
            dict: Departures[(i,j)] = {a: departure_time}
        """
        Departures = {}
        n_intervals = planning_horizon // interval
        transfer_nodes = sorted(C_minor)

        # === Forward direction (left → right) ===
        left_terminal = transfer_nodes[0]
        for j in range(1, len(transfer_nodes)):
            Departures[(left_terminal, transfer_nodes[j])] = {
                a: 150 + interval * a for a in range(int(n_intervals))
            }

        # Forward arcs (including skips)
        for i in range(1, len(transfer_nodes)):
            for j in range(i+1, len(transfer_nodes)):
                Departures[(transfer_nodes[i], transfer_nodes[j])] = {
                    a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k+1]] for k in range(i, j)))
                    for a, d in Departures[(left_terminal, transfer_nodes[i])].items()
                }

        # === Backward direction (right → left) ===
        right_terminal = transfer_nodes[-1]
        for j in range(len(transfer_nodes)-1):
            Departures[(right_terminal, transfer_nodes[j])] = {
                a: 150 + interval * a for a in range(int(n_intervals))
            }

        # Backward arcs (including skips)
        for i in range(len(transfer_nodes)-2, -1, -1):
            for j in range(i-1, -1, -1):
                Departures[(transfer_nodes[i], transfer_nodes[j])] = {
                    a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k-1]] for k in range(i, j, -1)))
                    for a, d in Departures[(right_terminal, transfer_nodes[i])].items()
                }

        return Departures

    ##### Arc Elimination #####

    def build_tij(self, t_transfer, nodes, n_requests, n_vehicles, n_trans_nodes,
                include_charging=True, duplicate=True, eliminate_arcs=False,
                di=None, ei=None, li=None, Lbar=None, P=None, D=None, N=None,
                zeroDepot=None, endDepot=None, Cr=None, n=None, pair_pi_di=None, pair_pi_di_M=None):
        
        """
        Build tij and cij with options for charging, duplication, and arc elimination.

        Args:
            t (np.array): base travel time matrix
            nodes (np.array): nodes array
            n_requests (int): number of requests
            n_vehicles (int): number of vehicles
            n_trans_nodes (int): number of physical transfer nodes
            include_charging (bool): include charging stations
            duplicate (bool): duplicate transfer and charging nodes per request/vehicle
            eliminate_arcs (bool): apply arc elimination rules
            di, ei, li, Lbar, P, D, N, zeroDepot, endDepot, Cr, n: optional sets/params 
                needed for arc elimination.

        Returns:
            tij (dict): travel time dictionary
            cij (dict): cost dictionary
        """

        # === Step 2. Build tij, cij dicts ===
        tij = {(self.base(i), self.base(j)): float(t_transfer[self.base(i), self.base(j)]) for i in N for j in N}

        # === Step 3. Arc elimination (optional) ===
        if eliminate_arcs and all(v is not None for v in [di, ei, li, Lbar, P, D, N, zeroDepot, endDepot, n]):

            # Depot restrictions
            for i in N:
                if (i, zeroDepot) in tij:
                    del tij[(self.base(i), self.base(zeroDepot))]
                if (endDepot, i) in tij:
                    del tij[(self.base(endDepot), self.base(i))]

            for i in D:
                if (zeroDepot, i) in tij:
                    del tij[(self.base(zeroDepot), self.base(i))]
            for i in P:
                if (i, endDepot) in tij:
                    del tij[(self.base(i), self.base(endDepot))]

            # Remove loops
            for i in N:
                if (i, i) in tij:
                    del tij[(self.base(i), self.base(i))]

            # Drop-off → pickup arcs
            for p, d in pair_pi_di.items():
                if (d, p) in tij:
                    del tij[(self.base(d), self.base(p))]

            for p, d in pair_pi_di_M.items():
                if (d, p) in tij:
                    del tij[(self.base(d), self.base(p))]

            # Time-window infeasibility
            for i in N:
                for j in N:
                    if (i, j) in tij and i in ei and i in di and j in li:
                        if ei[self.base(i)] + di[self.base(i)] + tij[(self.base(i), self.base(j))] >= li[self.base(j)]:
                            del tij[(self.base(i), self.base(j))]

            # Ride time infeasibility
            for i in P:
                for j in N:
                    if (self.base(i), self.base(j)) in tij and (self.base(j), self.base(pair_pi_di[i])) in tij:
                        if tij[(self.base(i), self.base(j))] + di[self.base(j)] + tij[(self.base(j), self.base(pair_pi_di[i]))] >= Lbar[self.base(i)]:
                            del tij[(self.base(i), self.base(j))]
                            del tij[(self.base(j), self.base(pair_pi_di[i]))]

            # Path infeasibility
            for i in P:
                for j in P:
                    if i != j:
                        if (self.base(j), self.base(i)) in tij and (self.base(i), self.base(pair_pi_di[j])) in tij and (self.base(pair_pi_di[j]), self.base(pair_pi_di[i])) in tij:
                            if tij[(self.base(j), self.base(i))] + di[self.base(i)] + tij[(self.base(i), self.base(pair_pi_di[i]))] + di[self.base(pair_pi_di[i])] + tij[(self.base(pair_pi_di[j]), self.base(pair_pi_di[i]))] > Lbar[self.base(i)]:
                                del tij[(self.base(i), self.base(pair_pi_di[i]))]
            for i in P:
                for j in P:
                    if i != j:
                        if (self.base(i), self.base(pair_pi_di[i])) in tij and (self.base(pair_pi_di[i]), self.base(j)) in tij and (self.base(j), self.base(pair_pi_di[i])) in tij:
                            if tij[(self.base(i), self.base(pair_pi_di[i]))] + di[self.base(pair_pi_di[i])] + tij[(self.base(pair_pi_di[i]), self.base(j))] + di[self.base(j)] + tij[(self.base(j), self.base(pair_pi_di[i]))] > Lbar[self.base(i)]:
                                del tij[(self.base(pair_pi_di[i]), self.base(j))]

            # Transfer-specific cleaning only if Cr exists
            if Cr is not None:
                for i in P:
                    drop_node = pair_pi_di[i]
                    if i in Cr:
                        for Si in Cr[i]:
                            if (drop_node, Si) in tij:
                                del tij[(self.base(drop_node), self.base(Si))]
                            if (Si, i) in tij:
                                del tij[(self.base(Si), self.base(i))]

        cij = tij.copy()
        return tij, cij

    ### === USAGE ===

    def pack_data(self):

        # Parametres

        ### Base travel time matrix
        if self. MoPS:
            # Node order:
            # 0 depot0, 1 P1, 2 P2, 3 P3, 4 P4, 5 D1, 6 D2, 7 D3, 8 D4,
            # 9 MoPS_pickup, 10 MoPS_dropoff, 11 depot1, 12 T1, 13 T2, 14 T3

            t = np.array([
                [  0,  23, 173, 131, 195, 170, 129, 160, 153, 157,  80,   0, 120,  21, 163],
                [ 23,   0, 166, 123, 173, 179, 125, 138, 168, 159, 103,  23, 135,  26, 159],
                [173, 166,   0,  57, 289,  41, 268, 254,  79,  59, 253, 173,  57, 191, 301],
                [131, 123,  57,   0, 258,  72, 236, 222,  60,  63, 211, 131,  29, 148, 270],
                [195, 173, 289, 258,   0, 318,  73,  35, 306, 293, 115, 195, 273, 174,  39],
                [170, 179,  41,  72, 318,   0, 296, 282,  53,  55, 250, 170,  51, 191, 329],
                [129, 125, 268, 236,  73, 296,   0,  48, 282, 270, 128, 129, 248, 108,  34],
                [160, 138, 254, 222,  35, 282,  48,   0, 271, 257,  80, 160, 237, 139,  57],
                [153, 168,  79,  60, 306,  53, 282, 271,   0,  64, 233, 153,  35, 174, 316],
                [157, 159,  59,  63, 293,  55, 270, 257,  64,   0, 237, 157,  43, 176, 304],
                [ 80, 103, 253, 211, 115, 250, 128,  80, 233, 237,   0,  80, 200, 101, 137],
                [  0,  23, 173, 131, 195, 170, 129, 160, 153, 157,  80,   0, 120,  21, 163],
                [120, 135,  57,  29, 273,  51, 248, 237,  35,  43, 200, 120,   0, 141, 283],
                [ 21,  26, 191, 148, 174, 191, 108, 139, 174, 176, 101,  21, 141,   0, 142],
                [163, 159, 301, 270,  39, 329,  34,  57, 316, 304, 137, 163, 283, 142,   0],
            ], dtype=float)

            nodes = np.array([
                [0, "depot0", "", 0, 86400, 0, 0],
                [11, "depot1", "", 0, 86400, 0, 0],
                [1, "pickup", 1, 50, 950, 5, 1],
                [5, "dropoff", 1, "", "", 5, 0],
                [2, "pickup", 2, 350, 1250, 5, 1],
                [6, "dropoff", 2, "", "", 5, 0],
                [3, "pickup", 3, 650, 1550, 5, 1],
                [7, "dropoff", 3, "", "", 5, 0],
                [4, "pickup", 4, 950, 1850, 5, 1],
                [8, "dropoff", 4, "", "", 5, 0],
                [9, "MoPS pickup", 5, 500, 1400, 5, 1],
                [10, "MoPS dropoff", 5, "", "", 5, 0],
            ], dtype=object)

        else:
            t = np.array([
                [  0,  23, 173, 131, 195, 170, 129, 160, 153,   0, 120,  21, 163],
                [ 23,   0, 166, 123, 173, 179, 125, 138, 168,  23, 135,  26, 159],
                [173, 166,   0,  57, 289,  41, 268, 254,  79, 173,  57, 191, 301],
                [131, 123,  57,   0, 258,  72, 236, 222,  60, 131,  29, 148, 270],
                [195, 173, 289, 258,   0, 318,  73,  35, 306, 195, 273, 174,  39],
                [170, 179,  41,  72, 318,   0, 296, 282,  53, 170,  51, 191, 329],
                [129, 125, 268, 236,  73, 296,   0,  48, 282, 129, 248, 108,  34],
                [160, 138, 254, 222,  35, 282,  48,   0, 271, 160, 237, 139,  57],
                [153, 168,  79,  60, 306,  53, 282, 271,   0, 153,  35, 174, 316],
                [  0,  23, 173, 131, 195, 170, 129, 160, 153,   0, 120,  21, 163],
                [120, 135,  57,  29, 273,  51, 248, 237,  35, 120,   0, 141, 283],
                [ 21,  26, 191, 148, 174, 191, 108, 139, 174,  21, 141,   0, 142],
                [163, 159, 301, 270,  39, 329,  34,  57, 316, 163, 283, 142,   0]
            ], dtype=float)

            nodes = np.array([
                [0, "depot0", "", 0, 86400, 0, 0],
                [9, "depot1", "", 0, 86400, 0, 0],
                [1, "pickup", 1, 50, 950, 5, 1],
                [5, "dropoff", 1, "", "", 5, 0],
                [2, "pickup", 2, 350, 1250, 5, 1],
                [6, "dropoff", 2, "", "", 5, 0],
                [3, "pickup", 3, 650, 1550, 5, 1],
                [7, "dropoff", 3, "", "", 5, 0],
                [4, "pickup", 4, 950, 1850, 5, 1],
                [8, "dropoff", 4, "", "", 5, 0],
            ], dtype=object)

        # print("nodes:", nodes)

        n_requests = len(set(nodes[nodes[:,2] != "", 2].astype(int)))
        n_vehicles = 2 
        n_trans_nodes = len(t) - len(nodes)

        # print("n_requests:", n_requests)
        # print("n_trans_nodes:", n_trans_nodes)

        Departures = None
        Cr = None

        # Build matrices and nodes
        t_transfer = self.build_t_transfer(t, n_requests, n_vehicles, n_trans_nodes, nodes)
        # print("t_transfer: ",t_transfer)
        nodes = self.build_nodes(nodes, n_requests, n_vehicles, n_trans_nodes)
        # print("nodes after nuild_nodes :", nodes)

        # Build sets Id
        # N, P, D, C, F, zeroDepot, endDepot, N_type = build_node_id_sets(nodes)
        if self.use_imjn:
            R = self.build_requests(nodes)
            N, P, D, C, F, zeroDepot, endDepot, P_M, D_M= self.build_node_imjn_sets(nodes)
            N_minor, P_minor, D_minor, C_minor, F_minor, zeroDepot_bin, endDepot_bin, P_M_minor, D_M_minor= self.build_node_id_sets(nodes)
            Departures = self.build_departures(t_transfer, C_minor)
        else: 
            R, Cr = self.build_requests(nodes)
            N, P, D, C, F, zeroDepot, endDepot, P_M, D_M = self.build_node_id_sets(nodes)

        fi_r = self.build_fi_r(nodes, N, R, P, D, P_M, D_M)
        # print("fi_r is: ",fi_r)

        # Build parametres
        di, ei, li, qr, pair_pi_di, Lbar, pair_pi_di_M = self.build_service_params(nodes, t_transfer, P, D, P_M, D_M)

        # Constants
        Q = 4
        T = 3600 * 4
        V = n_vehicles
        K = list(range(V))
        ek = {k: 0 for k in K}
        lk = {k: 86400 for k in K}
        M = 100000
        n = len(P)
        # Electric Vehicle Constants
        C_bat_kWh = 100 # Battery Capacity in kWh
        C_bat = C_bat_kWh*(60*60) # Battery Capacity in kWs
        alpha = 200 # Charging Power in kW
        beta = 15 # Average Power in kW
        gamma = 0.1 # Minimum battery capacity set at 10%

        # --- Dataset with all variants included ---

        sets = dict(
            # === Node & Set Information ===
            nodes=nodes, N=N, P=P, D=D, C=C, F=F, R=R, K=K,
            zeroDepot=zeroDepot, endDepot=endDepot, P_M=P_M, D_M=D_M)
        params = dict(
            # === Parameters ===
            Q=Q, T=T, qr=qr, di=di, ei=ei, li=li,
            ek=ek, lk=lk, fi_r=fi_r, Lbar=Lbar, pair_pi_di=pair_pi_di, pair_pi_di_M=pair_pi_di_M, M=M, n=n)

        params['tij'], params['cij'] = self.build_tij(
            t_transfer, nodes, n_requests, n_vehicles, n_trans_nodes,
            include_charging=self.ev_constraints, duplicate=self.duplicate_transfers,
            eliminate_arcs=self.arc_elimination,
            di=di, ei=ei, li=li, Lbar=Lbar, P=P, D=D, N=N,
            zeroDepot=zeroDepot, endDepot=endDepot, Cr=Cr if self.duplicate_transfers else None,
            n=n, pair_pi_di=pair_pi_di, pair_pi_di_M=pair_pi_di_M)
        
        tij_keys = params['tij'].keys()
        sets["A"] = {(i, j) for i in N for j in N if (self.base(i), self.base(j)) in tij_keys and self.base(i) != self.base(j)}

        if self.duplicate_transfers:
            params["Cr"] = Cr 
        else: params["Cr"] = None
        if self.ev_constraints:
            params["C_bat"] = C_bat
            params["alpha"] = alpha
            params["beta"] = beta
            params["gamma"] = gamma
        if self.use_imjn:
            params['Departures'] = Departures
        else: params['Departures'] = None

        return sets, params


    def build(self):
        self.sets, self.params = self.pack_data()
        return self.sets, self.params
    
if __name__ == "__main__":
    builder = DARPDataBuilder(
        duplicate_transfers=True,
        arc_elimination=True,
        ev_constraints=False,
        use_imjn=False,
        MoPS=True
    )
    sets, params = builder.build()
    print("✅ Build complete.")
    print("Node definition:", sets['nodes'])
    print("Number of nodes:", len(sets["N"]))
    print("Number of arcs:", len(sets["A"]))
    print("Travel time dict size:", len(params["tij"]))
    # print("f_ir:", params['fi_r'])
    # print("tij: ", params['tij'])
    print("ei:", params['ei'])
    print("li:", params['li'])
    print("Example C:", list(sets['C'])[:5])


