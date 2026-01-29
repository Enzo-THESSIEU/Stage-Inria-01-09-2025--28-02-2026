##### Imports #####
import numpy as np
from extract_data import data_extractor
from sql_to_instance import build_random_instance

class DARPDataBuilder:
    def __init__(self, duplicate_transfers=True, arc_elimination=True,
                 ev_constraints=False, use_imjn=False, MoPS = False, datafile_instance=False, max_visits_transfer=5):
        self.duplicate_transfers = duplicate_transfers
        self.arc_elimination = arc_elimination
        self.ev_constraints = ev_constraints
        self.use_imjn = use_imjn
        self.MoPS = MoPS
        self.datafile_instance = datafile_instance
        self.max_visits_transfer = max_visits_transfer

    def base(self, x):
        return x[0] if isinstance(x, tuple) else x
        
    def build_t_transfer(self, t, n_requests, n_vehicles, n_trans_nodes):
        """
        Build the extended travel-time matrix with:
        - transfer nodes (shared or duplicated per request)
        - charging nodes (duplicated per vehicle)
        - zero-distance between paired transfer/charging nodes

        Parameters
        ----------
        t : np.ndarray
            Base travel-time matrix (n_base × n_base)
        n_requests : int
        n_vehicles : int
        n_trans_nodes : int   # number of transfer nodes per cluster

        Returns
        -------
        t_ext : np.ndarray
            Extended matrix with base, transfers, and charging nodes
        """

        # === 1. BASIC SIZES ===
        n_base = t.shape[0]

        # Transfer nodes
        if self.duplicate_transfers:
            n_transfer_layers = n_requests - 1       # one layer per request
        else:
            n_transfer_layers = 1                   # one shared layer

        n_transfers = n_transfer_layers * n_trans_nodes

        # Charging nodes
        if self.ev_constraints:
            n_charging = n_vehicles * n_trans_nodes
        else:
            n_charging = 0

        # Total extended matrix size
        n_total = n_base + n_transfers + n_charging

        # Init extension
        t_ext = np.full((n_total, n_total), np.inf)
        np.fill_diagonal(t_ext, 0)

        # Copy base matrix
        t_ext[:n_base, :n_base] = t


        # === 2. INDEXING ===
        trans_start   = n_base
        trans_end     = trans_start + n_transfers

        charge_start  = trans_end
        charge_end    = charge_start + n_charging


        # === 3. TRANSFER NODE BLOCK ===
        if n_transfers > 0:

            # We duplicate the base transfer block t[T,T]
            # Extract the transfer block from the base matrix
            # Transfer nodes are assumed contiguous at the end of base structure
            t_T = t[n_base - n_trans_nodes : n_base, n_base - n_trans_nodes : n_base]

            # --- Transfer → Transfer (tiled by number of layers)
            t_ext[trans_start:trans_end, trans_start:trans_end] = \
                np.tile(t_T, (n_transfer_layers, n_transfer_layers))

            # --- Base → Transfer (tiled horizontally)
            t_base_to_T = t[:, n_base - n_trans_nodes : n_base]
            t_ext[:n_base, trans_start:trans_end] = \
                np.tile(t_base_to_T, n_transfer_layers)

            # --- Transfer → Base (tiled vertically)
            t_T_to_base = t[n_base - n_trans_nodes : n_base, :]
            t_ext[trans_start:trans_end, :n_base] = \
                np.tile(t_T_to_base, (n_transfer_layers, 1))


        # === 4. CHARGING NODES BLOCK ===
        if n_charging > 0:

            # Charging nodes = transfer nodes duplicated per vehicle
            # Each transfer node k gives charging nodes: k1, k2, ..., kV

            # --- Transfer ↔ Charging: zero distance (same location)
            for v in range(n_vehicles):
                offset = charge_start + v * n_trans_nodes

                # zero distance transfer_i <-> charging_{i,v}
                for i in range(n_transfers):
                    t_ext[trans_start + i, offset + (i % n_trans_nodes)] = 0
                    t_ext[offset + (i % n_trans_nodes), trans_start + i] = 0

            # --- Charging ↔ Charging block (copy transfer distances)
            t_ext[charge_start:charge_end, charge_start:charge_end] = \
                np.tile(t_T, (n_vehicles, n_vehicles))

            # --- Base → Charging
            t_base_to_T = t[:, n_base - n_trans_nodes : n_base]
            t_ext[:n_base, charge_start:charge_end] = \
                np.tile(t_base_to_T, n_vehicles)

            # --- Charging → Base
            t_T_to_base = t[n_base - n_trans_nodes : n_base, :]
            t_ext[charge_start:charge_end, :n_base] = \
                np.tile(t_T_to_base, (n_vehicles, 1))


        return t_ext
    
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

        if self.datafile_instance:
            return base_nodes

        transfers, charging_stations = [], []

        # --- TRANSFER NODES ---
        transfer_start_id = max(set(base_nodes[:, 0])) + 1
        if self.duplicate_transfers:
            # One set of transfers per request
            for r in range(1, n_requests + 1):
                for t_id in range(n_trans_nodes):
                    node_id = transfer_start_id + (r - 1) * n_trans_nodes + t_id
                    transfers.append([node_id, "transfer", str(r), 0, 86400, 5, 0])
        else:
            # Single shared set of transfers
            for t_id in range(n_trans_nodes):
                node_id = transfer_start_id + t_id
                transfers.append([node_id, "transfer", "", 0, 86400, 5, 0])

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
        # print("N", N)
        P = nodes[nodes[:, 1] == "pickup", 0].astype(int).tolist()
        D = nodes[nodes[:, 1] == "dropoff", 0].astype(int).tolist()
        C = nodes[nodes[:, 1] == "transfer", 0].astype(int).tolist()
        if self.datafile_instance:
            zeroDepot = nodes[nodes[:, 1] == "depot0", 0].astype(int).tolist()
            endDepot = nodes[nodes[:, 1] == "depot1", 0].astype(int).tolist()
        else:
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
        if self.datafile_instance:
            zeroDepot_node = [(zeroDepot, 1) for zeroDepot in zeroDepot_node]
            endDepot_node = [(endDepot, 2) for endDepot in endDepot_node]

        else:
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
        if self.datafile_instance:
            N = zeroDepot_node + endDepot_node + P + D + C + F + P_M + D_M
        else:
            N = [zeroDepot_node, endDepot_node] + P + D + C + F + P_M + D_M

        # print("P", P)
        # print("D", D)

        return N, P, D, C, F, zeroDepot_node, endDepot_node, P_M, D_M

    def build_requests(self, nodes):
        if self.use_imjn:
            r_values = nodes[np.isin(nodes[:, 1], ["pickup", "dropoff", "transfer", "MoPS pickup", "MoPS dropoff"]), 2]
            r_values = r_values[r_values != ""]
            R = sorted(set(r_values.astype(int)))
            # print("RR", R)
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
            if not self.datafile_instance:
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
            # print("self.base(p)", self.base(p))
            # print("nodes[:, 0]", int(nodes[:, 0]))
            # print("nodes[nodes[:, 0]]", nodes[nodes[:, 0]])
            # print("nodes[nodes[:, 0] == self.base(p), 2", nodes[nodes[:, 0] == self.base(p), 2])
            if self.datafile_instance:
                r = int(nodes[nodes[:, 0].astype(int) == self.base(p), 2][0])
            else:
                r = int(nodes[nodes[:, 0] == self.base(p), 2][0])
            fi_r[r, p] = 1
            fi_r[r, d] = -1
        if self.MoPS:
            for p, d in zip(P_M, D_M):
                r = int(nodes[nodes[:, 0] == self.base(p), 2][0])
                fi_r[r, p] = 1
                fi_r[r, d] = -1
        return fi_r

    def build_departures(self, t, C_minor, interval=10, planning_horizon=720):    ##### Modify again ######
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
                a: interval * a for a in range(int(n_intervals))
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
                a: interval * a for a in range(int(n_intervals))
            }

        # Backward arcs (including skips)
        for i in range(len(transfer_nodes)-2, -1, -1):
            for j in range(i-1, -1, -1):
                Departures[(transfer_nodes[i], transfer_nodes[j])] = {
                    a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k-1]] for k in range(i, j, -1)))
                    for a, d in Departures[(right_terminal, transfer_nodes[i])].items()
                }

        # if self.MoPS:
        for i in C_minor:
            Departures[(i, i)] = {a: interval * a for a in range(int(n_intervals))}
            # Departures[(13,13)] = {a: interval * a for a in range(int(n_intervals))}
            # Departures[(14,14)] = {a: interval * a for a in range(int(n_intervals))}
            # Departures[(12,12)] = {a: interval * a for a in range(int(n_intervals))}

        # else:
        #     Departures[(10,10)] = {a: interval * a for a in range(int(n_intervals))}
        #     Departures[(11,11)] = {a: interval * a for a in range(int(n_intervals))}
        #     Departures[(12,12)] = {a: interval * a for a in range(int(n_intervals))}

        return Departures

    ##### Arc Elimination #####

    def build_tij(self, t_transfer, nodes, n_requests, n_vehicles, n_trans_nodes,
                di=None, ei=None, li=None, Lbar=None, P=None, D=None, N=None,
                zeroDepot=None, endDepot=None, Cr=None, n=None, pair_pi_di=None, pair_pi_di_M=None,
                C=None, R=None, P_M=None, D_M=None):
        
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
        # print(tij.keys())


        # === Step 3. Arc elimination (optional) ===
        if self.arc_elimination and all(v is not None for v in [di, ei, li, Lbar, P, D, N, zeroDepot, endDepot, n]):

            # Helper for short access
            base = self.base  

            # # --- Depot restrictions ---
            # if self.datafile_instance:
            #     for i in N:
            #         bi = base(i)
            #         for Depot_node in zeroDepot:   ### for this type of instance the zero depot node is the same as the end depot node
            #             if (bi, base(Depot_node)) in tij:
            #                 del tij[(bi, base(Depot_node))]
            #             if (base(Depot_node), bi) in tij:
            #                 del tij[(base(Depot_node), bi)]

            #     for i in D:
            #         for Depot_node in zeroDepot:
            #             if (base(Depot_node), base(i)) in tij:
            #                 del tij[(base(Depot_node), base(i))]
            #     for i in P:
            #         for Depot_node in zeroDepot:
            #             if (base(i), base(Depot_node)) in tij:
            #                 del tij[(base(i), base(Depot_node))]
                
            #     if self.MoPS:
            #         for i in D_M:
            #             for Depot_node in zeroDepot:
            #                 if (base(Depot_node), base(i)) in tij:
            #                     del tij[(base(Depot_node), base(i))]
            #         for i in P_M:
            #             for Depot_node in zeroDepot:
            #                 if (base(i), base(Depot_node)) in tij:
            #                     del tij[(base(i), base(Depot_node))]

            # else:
            for i in N:
                bi = base(i)
                # print("bi:", bi, type(bi))
                # print("zeroDepot:", zeroDepot, type(zeroDepot))
                # print("base(zeroDepot):", base(zeroDepot), type(base(zeroDepot)))
                if (bi, base(zeroDepot)) in tij:
                    del tij[(bi, base(zeroDepot))]
                if (base(endDepot), bi) in tij:
                    del tij[(base(endDepot), bi)]

            for i in D:
                if (base(zeroDepot), base(i)) in tij:
                    del tij[(base(zeroDepot), base(i))]
            for i in P:
                if (base(i), base(endDepot)) in tij:
                    del tij[(base(i), base(endDepot))]
            
            if self.MoPS:
                for i in D_M:
                    if (base(zeroDepot), base(i)) in tij:
                        del tij[(base(zeroDepot), base(i))]
                for i in P_M:
                    if (base(i), base(endDepot)) in tij:
                        del tij[(base(i), base(endDepot))]

            # # --- Remove self-loops ---
            # for i in N:
            #     if (base(i), base(i)) in tij:
            #         del tij[(base(i), base(i))]

            # --- Drop-off → Pickup arcs ---
            for p, d in pair_pi_di.items():
                if (base(d), base(p)) in tij:
                    del tij[(base(d), base(p))]

            if self.MoPS:
                for p, d in pair_pi_di_M.items():
                    if (base(d), base(p)) in tij:
                        del tij[(base(d), base(p))]

            # --- Time-window infeasibility ---
            # for i in N:
            #     for j in N:
            #         bi, bj = base(i), base(j)
            #         if (bi, bj) in tij and bi in ei and bi in di and bj in li:
            #             if ei[bi] + di[bi] + tij[(bi, bj)] >= li[bj]:
            #                 del tij[(bi, bj)]

            # --- Ride time infeasibility ---
            # for i in P:
            #     for j in N:
            #         bi, bj = base(i), base(j)
            #         if (bi, bj) in tij and (bj, base(pair_pi_di[i])) in tij:
            #             if tij[(bi, bj)] + di[bj] + tij[(bj, base(pair_pi_di[i]))] >= Lbar[bi]:
            #                 if (bi, bj) in tij:
            #                     del tij[(bi, bj)]
            #                 if (bj, base(pair_pi_di[i])) in tij:
            #                     del tij[(bj, base(pair_pi_di[i]))]

            if self.MoPS:
                for i in P_M:
                    for j in N:
                        bi, bj = base(i), base(j)
                        if (bi, bj) in tij and (bj, base(pair_pi_di_M[i])) in tij:
                            if tij[(bi, bj)] + di[bj] + tij[(bj, base(pair_pi_di_M[i]))] >= Lbar[bi]:
                                if (bi, bj) in tij:
                                    del tij[(bi, bj)]
                                if (bj, base(pair_pi_di[i])) in tij:
                                    del tij[(bj, base(pair_pi_di[i]))]


            # --- Path infeasibility ---
            pair_all = {**pair_pi_di, **pair_pi_di_M}

            if self.MoPS:
                for i in P + P_M:
                    for j in P + P_M:
                        if i != j:
                            bi, bj = base(i), base(j)
                            dpi, dpj = base(pair_all[i]), base(pair_all[j])
                            if (bj, bi) in tij and (bi, dpi) in tij and (dpj, dpi) in tij:
                                if tij[(bj, bi)] + di[bi] + tij[(bi, dpi)] + di[dpi] + tij[(dpj, dpi)] > Lbar[bi]:
                                    if (bi, dpi) in tij:
                                        del tij[(bi, dpi)]
                for i in P + P_M:
                    for j in P + P_M:
                        if i != j:
                            bi, bj = base(i), base(j)
                            dpi = base(pair_all[i])
                            if (bi, dpi) in tij and (dpi, bj) in tij and (bj, dpi) in tij:
                                if tij[(bi, dpi)] + di[dpi] + tij[(dpi, bj)] + di[bj] + tij[(bj, dpi)] > Lbar[bi]:
                                    if (dpi, bj) in tij:
                                        del tij[(dpi, bj)]

            else:
                for i in P:
                    for j in P:
                        if i != j:
                            bi, bj = base(i), base(j)
                            dpi, dpj = base(pair_pi_di[i]), base(pair_pi_di[j])
                            if (bj, bi) in tij and (bi, dpi) in tij and (dpj, dpi) in tij:
                                if tij[(bj, bi)] + di[bi] + tij[(bi, dpi)] + di[dpi] + tij[(dpj, dpi)] > Lbar[bi]:
                                    if (bi, dpi) in tij:
                                        del tij[(bi, dpi)]
                for i in P:
                    for j in P:
                        if i != j:
                            bi, bj = base(i), base(j)
                            dpi = base(pair_pi_di[i])
                            if (bi, dpi) in tij and (dpi, bj) in tij and (bj, dpi) in tij:
                                if tij[(bi, dpi)] + di[dpi] + tij[(dpi, bj)] + di[bj] + tij[(bj, dpi)] > Lbar[bi]:
                                    if (dpi, bj) in tij:
                                        del tij[(dpi, bj)]



            # --- Remove dropoff→transfer and transfer→pickup arcs (per request) ---
            if Cr is not None:
                if self.MoPS:
                    for p in P + P_M:
                        drop_node = pair_all[p]
                        if base(p) in Cr:   # safer lookup by base
                            for Si in Cr[base(p)]:
                                if (base(drop_node), base(Si)) in tij:
                                    del tij[(base(drop_node), base(Si))]
                                if (base(Si), base(p)) in tij:
                                    del tij[(base(Si), base(p))]

                else:
                    for p in P:
                        drop_node = pair_pi_di[p]
                        if base(p) in Cr:   # safer lookup by base
                            for Si in Cr[base(p)]:
                                if (base(drop_node), base(Si)) in tij:
                                    del tij[(base(drop_node), base(Si))]
                                if (base(Si), base(p)) in tij:
                                    del tij[(base(Si), base(p))]


            # # # Remove arcs from different transfer nodes of different requests
            # # if Cr is not None:
            # #     for r in R:
            # #         for i in C:
            # #             for j in Cr[r]:
            # #                 if i not in Cr[r]:
            # #                     if (self.base(i), self.base(j)) in tij:
            # #                         del tij[self.base(i), self.base(j)]
            # #                     if (self.base(j), self.base(i)) in tij:
            # #                         del tij[self.base(j), self.base(i)]


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
        t_transfer = self.build_t_transfer(t, n_requests, n_vehicles, n_trans_nodes)
        # print("t_transfer: ",t_transfer)
        nodes = self.build_nodes(nodes, n_requests, n_vehicles, n_trans_nodes)
        # print("nodes after nuild_nodes :", nodes)

        if self.datafile_instance:
            # file_path_prev = "C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes\\Stage-Inria-01-09-2025--28-02-2026\\code and instances\\"
            # extractor = data_extractor(file_path = file_path_prev + "datafile0.txt", mrt_factor = 2)
            # nodes = extractor.extract_node_info()
            # t = extractor.extract_travel_matrix(pt_speed = 1)
            # mrt = extractor.compute_maximum_ride_time(t)
            # nodes = np.array(extractor.tighten_time_windows(nodes, t))
            # n_vehicles = extractor.n_vehicles
            # n_requests = extractor.n_requests
            # n_transfers = extractor.n_transfers
            # t_transfer = t
            nodes, t = build_random_instance()
            nodes = np.array(nodes, dtype=object)
            t = np.array(t, dtype=float)
            n_vehicles = 4
            n_requests = 10
            t_transfer = t

        # Build sets Id
        # N, P, D, C, F, zeroDepot, endDepot, N_type = build_node_id_sets(nodes)
        if self.use_imjn:
            R = self.build_requests(nodes)
            N, P, D, C, F, zeroDepot, endDepot, P_M, D_M= self.build_node_imjn_sets(nodes)
            print("zero Depot", zeroDepot)
            print("P", P)
            N_minor, P_minor, D_minor, C_minor, F_minor, zeroDepot_bin, endDepot_bin, P_M_minor, D_M_minor= self.build_node_id_sets(nodes)
            Departures = self.build_departures(t_transfer, C_minor)
        else: 
            R, Cr = self.build_requests(nodes)
            N, P, D, C, F, zeroDepot, endDepot, P_M, D_M = self.build_node_id_sets(nodes)

        if self.datafile_instance:
            zeroDepot = zeroDepot[0]
            endDepot = endDepot[0]

        # print("R fi_r", R)
        # print("N fi_r", N)
        fi_r = self.build_fi_r(nodes, N, R, P, D, P_M, D_M)
        # print("fi_r is: ",fi_r)

        # Build parametres
        di, ei, li, qr, pair_pi_di, Lbar, pair_pi_di_M = self.build_service_params(nodes, t_transfer, P, D, P_M, D_M)
        print("ei :", ei)
        print("li :", li)
        # Constants
        Q = 4
        T = 3600 * 4
        V = n_vehicles
        K = list(range(V))
        ek = {k: 0 for k in K}
        # lk = {k: 86400 for k in K}
        lk = {k: 1440 for k in K}
        M = 100000
        n = len(P)
        # Electric Vehicle Constants
        C_bat_kWh = 14.85 # Battery Capacity in kWh 
        C_bat = C_bat_kWh*(60*60) # Battery Capacity in kWs
        alpha = 2.97 * 10 # Charging Power in kW (Recharge complète en 5 heures donc puissance de recharge de 14,85kWh/5 heures = 2,97 kW)
        beta = 0.055 * 60 * 10  # Average Power in kW (0,055 kWh par minute = 3,3 kW)
        gamma = 0.1 # Minimum battery capacity set at 10%
        gamma_end = 0.7 # Minimum battery Capacity at end Depot

        def pair_pi_di_base_build(pair_pi_di):
            base_pair_pi_di = {}
            for key, val in pair_pi_di.items():
                base_pair_pi_di[self.base(key)] = self.base(val)
            return base_pair_pi_di

        base_pair_pi_di = pair_pi_di_base_build(pair_pi_di)
        base_pair_pi_di_M = pair_pi_di_base_build(pair_pi_di_M)

        # --- Dataset with all variants included ---

        sets = dict(
            # === Node & Set Information ===
            nodes=nodes, N=N, P=P, D=D, C=C, F=F, R=R, K=K,
            zeroDepot=zeroDepot, endDepot=endDepot, P_M=P_M, D_M=D_M)
        params = dict(
            # === Parameters ===
            Q=Q, T=T, qr=qr, di=di, ei=ei, li=li,
            ek=ek, lk=lk, fi_r=fi_r, Lbar=Lbar, pair_pi_di=pair_pi_di, pair_pi_di_M=pair_pi_di_M, M=M, n=n, 
            base_pair_pi_di = base_pair_pi_di, base_pair_pi_di_M = base_pair_pi_di_M)

        if self.MoPS:
            params['tij'], params['cij'] = self.build_tij(
                t_transfer, nodes, n_requests, n_vehicles, n_trans_nodes,
                di=di, ei=ei, li=li, Lbar=Lbar, P=P, D=D, N=N,
                zeroDepot=zeroDepot, endDepot=endDepot, Cr=Cr if self.duplicate_transfers else None,
                n=n, pair_pi_di=pair_pi_di, pair_pi_di_M=pair_pi_di_M, C=C, R=R, P_M=P_M, D_M=D_M)
        else:
            params['tij'], params['cij'] = self.build_tij(
                t_transfer, nodes, n_requests, n_vehicles, n_trans_nodes,
                di=di, ei=ei, li=li, Lbar=Lbar, P=P, D=D, N=N,
                zeroDepot=zeroDepot, endDepot=endDepot, Cr=Cr if self.duplicate_transfers else None,
                n=n, pair_pi_di=pair_pi_di, pair_pi_di_M=pair_pi_di_M, C=C, R=R)
            
        tij_keys = params['tij'].keys()
        sets["A"] = {(i, j) for i in N for j in N if (self.base(i), self.base(j)) in tij_keys and i != j}

        if self.duplicate_transfers:
            params["Cr"] = Cr 
        else: params["Cr"] = None
        if self.ev_constraints:
            params["C_bat"] = C_bat
            params["alpha"] = alpha
            params["beta"] = beta
            params["gamma"] = gamma
            params['gamma_end'] = gamma_end
        if self.use_imjn:
            params['Departures'] = Departures
        else: params['Departures'] = None

        return sets, params


    def build(self):
        self.sets, self.params = self.pack_data()
        return self.sets, self.params
    
if __name__ == "__main__":
    builder = DARPDataBuilder(
        duplicate_transfers=False,
        arc_elimination=True,
        ev_constraints=False,
        use_imjn=True,
        MoPS=False,
        datafile_instance=True
    )
    sets, params = builder.build()
    print("✅ Build complete.")
    print("Node definition:")
    print("ei :", params['ei'])
    print("li :", params['li'])
    # print(params['tij'].keys())
    # for node in sets['nodes']:
    #     print("====================\n")
    #     print(f"Node {node} \n")

    # print("Number of nodes:", len(sets["N"]))
    # print("Number of arcs:", len(sets["A"]))
    # print("Travel time dict:", params["tij"])
    # print("f_ir:", params['fi_r'])
    # print("tij: ", params['tij'])
    # print("ei:", params['ei'])
    # print("li:", params['li'])
    # print("C", sets['C'])
    # print("R", sets['R'])
    # def base(x):
    #     return x[0] if isinstance(x, tuple) else x
    # for (i,j) in sets['A']:
    #     if base(i) == base(j):
    #         print((i,j))
    # tij_keys = params['tij'].keys()
    # for key in tij_keys:
    #     if key[0] == 0:
    #         print(key)
    # for arc in sets['A']:
    #     if arc[0] == (0,1):
    #         print(arc)
    # print("N", sets['N'])
    # print("Departures", params['Departures'])

    # transfer_arcs = {(i,j) for (i,j) in sets['A'] if (i in sets['C'] and j in sets['C'])}
    # DAR_arcs = {(i,j) for (i,j) in sets['A'] if (i,j) not in transfer_arcs}
    # request_arcs = {(i,j) for (i,j) in DAR_arcs if not (i == sets['zeroDepot'] or j == sets['endDepot'])}
    # for (i, j) in request_arcs:
    #     t = 1
    #     print(i, j)
    # DAR_depot_arcs = {(i, j) for (i, j) in DAR_arcs if (i, j) not in request_arcs}
    # for (i, j) in transfer_arcs:
    #     print(i, j)
    print(params['base_pair_pi_di'])

