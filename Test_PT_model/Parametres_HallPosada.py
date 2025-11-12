##### Imports #####
import numpy as np

class DARPDataBuilder_HallPosada:
    """
    Build instance data for Model 2 (Häll & Posada 2017)
    – homogeneous fleet, single resource.
    Compatible with DARPConstraintBuilder_HallPosada.
    """

    def __init__(self):
        self.max_visits_transfer = 3   # maximum visits per transfer location

    def base(self, x):
        """Return base node id from (i,m) tuple or int."""
        return x[0] if isinstance(x, tuple) else x

    # ----------------------------------------------------------------------
    # Node + network structure
    # ----------------------------------------------------------------------
    def build_nodes(self):
        """Returns the base pickup, dropoff, depot, and transfer nodes."""
        nodes = np.array([
            [0, "depot0", "", 0, 86400, 5, 0],
            [9, "depot1", "", 0, 86400, 5, 0],
            [1, "pickup", 1, 50, 950, 5, 1],
            [5, "dropoff", 1, "", "", 5, 0],
            [2, "pickup", 2, 350, 1250, 5, 1],
            [6, "dropoff", 2, "", "", 5, 0],
            [3, "pickup", 3, 650, 1550, 5, 1],
            [7, "dropoff", 3, "", "", 5, 0],
            [4, "pickup", 4, 950, 1850, 5, 1],
            [8, "dropoff", 4, "", "", 5, 0],
        ], dtype=object)
        return nodes

    def build_travel_times(self):
        """Travel-time matrix identical to 4-request base dataset."""
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
        return t

    # ----------------------------------------------------------------------
    # Sets and parameters
    # ----------------------------------------------------------------------
    def build_sets(self, nodes, n_transfers=3, n_vehicles=2):
        """Build model sets: N, NP, ND, NG, K, R."""
        P = [int(x[0]) for x in nodes if x[1] == "pickup"]
        D = [int(x[0]) for x in nodes if x[1] == "dropoff"]
        R = sorted(set(int(r) for r in nodes[nodes[:,2] != "",2]))

        zeroDepot = int(nodes[nodes[:,1]=="depot0",0][0])
        endDepot  = int(nodes[nodes[:,1]=="depot1",0][0])

        transfer_ids = [10, 11, 12]
        NG = [(i, m) for i in transfer_ids for m in range(1, self.max_visits_transfer+1)]
        NP = [(i,1) for i in P]
        ND = [(i,1) for i in D]
        N  = [(zeroDepot,1), (endDepot,1)] + NP + ND + NG
        K = list(range(n_vehicles))

        return dict(N=N, NP=NP, ND=ND, NG=NG, K=K, R=R,
                    zeroDepot=(zeroDepot,1), endDepot=(endDepot,1))

    def build_departures(self, t, transfer_nodes, interval=300, horizon=3600*4):
        """Build departures dictionary between transfer nodes (forward+backward)."""
        Departures = {}
        n_intervals = int(horizon/interval)
        for (i,j) in [(10,11),(11,12),(12,11),(11,10)]:
            Departures[(i,j)] = {d: 650 + interval*d for d in range(n_intervals)}
        return Departures

    def build_service_params(self, nodes, t, NP, ND):
        """Basic service-time, load, and time-window parameters."""
        di = {int(r[0]): float(r[5]) if r[5] != "" else 5 for r in nodes}
        ei = {int(r[0]): float(r[3]) for r in nodes if r[3] != ""}
        li = {int(r[0]): float(r[4]) for r in nodes if r[4] != ""}
        pair_pi_di = {p[0]: d[0] for p, d in zip(NP, ND)}  # flat link
        Lr = {i[0]: 1 for i in NP}
        Q = 4
        M = 100000
        return di, ei, li, pair_pi_di, Lr, Q, M

    def build_costs_times(self, t):
        """Build dictionaries for Cij and Tij (same numeric values)."""
        indices = range(t.shape[0])
        Tij = {(i,j): float(t[i,j]) for i in indices for j in indices}
        Cij = Tij.copy()
        Cij_dr = {(i,j,d,r): 0.5*Tij[(i,j)] for i,j in [(10,11),(11,12),(12,11),(11,10)]
                  for d in range(3) for r in range(1,5)}
        return Cij, Cij_dr, Tij

    # ----------------------------------------------------------------------
    # Build complete instance
    # ----------------------------------------------------------------------
    def pack_data(self):
        nodes = self.build_nodes()
        t = self.build_travel_times()

        sets = self.build_sets(nodes)
        NP, ND, NG = sets["NP"], sets["ND"], sets["NG"]
        N, R, K = sets["N"], sets["R"], sets["K"]
        zeroDepot, endDepot = sets["zeroDepot"], sets["endDepot"]

        Departures = self.build_departures(t, [10,11,12])
        di, ei, li, pair_pi_di, Lr, Q, M = self.build_service_params(nodes, t, NP, ND)
        Cij, Cij_dr, Tij = self.build_costs_times(t)

        # Build arc list
        A = [(i,j) for i in N for j in N if i != j]

        # Request balance (F_ir)
        Fir = {}
        for r in R:
            for (i, m) in NP:
                Fir[(r, (i, m))] = 1 if int(r) == int(i) else 0
            for (i, m) in ND:
                # find matching pickup (inverse lookup)
                match_pickups = [p for p,d in pair_pi_di.items() if d == i]
                Fir[(r, (i, m))] = -1 if match_pickups and int(r) == match_pickups[0] else 0
            for (i, m) in NG:
                Fir[(r, (i, m))] = 0

        # ------------------------------------------------------------------
        # Parameter dictionary
        # ------------------------------------------------------------------
        params = dict(
            Cij=Cij,
            Cij_dr=Cij_dr,
            Tij=Tij,
            TD={(i,j,d): Departures[i,j][d] for (i,j) in Departures for d in Departures[i,j]},
            TA={(i,j,d): Departures[i,j][d] + Tij[(i,j)] for (i,j) in Departures for d in Departures[i,j]},
            Dij=Departures,
            di=di, ei=ei, li=li,
            pair_pi_di=pair_pi_di,
            Lr=Lr, Q=Q, M=M,
            Fir=Fir
        )

        # Add remaining mandatory parameters expected by constraints
        T0j = {j: 0 for (i,j) in A if self.base(i) == sets["zeroDepot"][0]}
        params["T0j"] = T0j
        params["Ti_e"] = ei
        params["Ti_l"] = li
        params["Bi"]   = {i: 1800 for i in [p[0] for p in NP]}  # ride time limit
        params["rbar"] = 4  # offset between pickup and dropoff indices

        sets["A"] = A
        return sets, params

    def build(self):
        return self.pack_data()

# ----------------------------------------------------------------------
# Execute test
# ----------------------------------------------------------------------
if __name__ == "__main__":
    builder = DARPDataBuilder_HallPosada()
    sets, params = builder.build()
    print("✅ Hall & Posada dataset built.")
    print("Nodes:", len(sets["N"]))
    print("Transfers:", len(sets["NG"]))
    print("Requests:", sets["R"])
    print("Arcs:", len(sets["A"]))
    print("Param keys:", list(params.keys()))
