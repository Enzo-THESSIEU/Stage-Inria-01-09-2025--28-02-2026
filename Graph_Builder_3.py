##### Imports #####
import numpy as np

##### Load Nodes and Arcs #####
### Base Matrix
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

# --- Build expanded matrix (22×22) ---
n_base = 13        # original nodes
n_requests = 4
n_trans_nodes = 3        # 3 transfer stops × 3 physical nodes (10,11,12)
n_trans = n_trans_nodes * (n_requests - 1)
n_vehicles = 2
n_charge_nodes = n_trans_nodes
n_charge = n_vehicles * n_charge_nodes       # 3 Charging facilities loacated at transfer nodes for 2 vehciles
n_total = n_base + n_trans + n_charge

t_transfer = np.zeros((n_total, n_total))

# 1. Copy base block
t_transfer[:n_base, :n_base] = t

# 2. Map base nodes to transfer nodes (rows)
# replicate cols 10:13 across each transfer group
t_transfer[:n_base, n_base:] = np.tile(t[:, 10:13], (1, (n_requests - 1 + n_vehicles)))

# 3. Map transfer nodes to base nodes (cols)
# replicate rows 10:13 across each transfer group
t_transfer[n_base:, :n_base] = np.tile(t[10:13, :], ((n_requests - 1 + n_vehicles), 1))

# 4. Transfer-to-transfer block
# replicate 10:13×10:13 into a 9×9 block (3 groups × 3 groups)
t_transfer[n_base:, n_base:] = np.tile(t[10:13, 10:13], ((n_requests - 1 + n_vehicles), (n_requests - 1 + n_vehicles)))

# print(t_transfer)

nodes = np.array([
    # id,   type,     r,   e,    l,     d,   q
    [0,   "depot0",  "",   0, 86400,   0,   0],
    [9,   "depot1",  "",   0, 86400,   0,   0],
    [1,   "pickup",   1,  50,   950,   5,   1],
    [5,   "dropoff",  1,  "",   "",    5,   0],
    [2,   "pickup",   2, 350,  1250,   5,   1],
    [6,   "dropoff",  2,  "",   "",    5,   0],
    [3,   "pickup",   3, 650,  1550,   5,   1],
    [7,   "dropoff",  3,  "",   "",    5,   0],
    [4,   "pickup",   4, 950,  1850,   5,   1],
    [8,   "dropoff",  4,  "",   "",    5,   0],
], dtype=object)

# Parameters
n_requests = 4
n_vehicles = 2
n_trans_nodes = 3   # physical transfer nodes (10,11,12)

# Build transfer nodes: duplicate 10..12 for each request
transfers = []
next_id = 10
for r in range(1, n_requests + 1):
    for phys_id in [10, 11, 12]:
        transfers.append([next_id, "transfer", str(r), "", "", 5, 0])
        next_id += 1

# Build charging stations: duplicate 10..12 for each vehicle
charging_stations = []
for k in range(1, n_vehicles + 1):
    for phys_id in [10, 11, 12]:
        charging_stations.append([next_id, "charging station", "", 0, 86400, 0, 0])
        next_id += 1

# Stack everything
nodes = np.vstack([
    nodes,
    np.array(transfers, dtype=object),
    np.array(charging_stations, dtype=object)
])

# print("Nodes:\n", nodes)

##### Build sets #####
# Nodes
N = nodes[:,0].astype(int).tolist()
N_type = nodes[:,1].astype(str).tolist()
P = nodes[nodes[:,1]=="pickup",0].astype(int).tolist()
D = nodes[nodes[:,1]=="dropoff",0].astype(int).tolist()
C = nodes[nodes[:,1]=="transfer",0].astype(int).tolist()
F = nodes[nodes[:,1]=="charging station",0].astype(int).tolist()
# D = [300 * i for i in range(24*60*60/300)] 
# print(C)
# print(F)



# print(len(N), len(P), len(D), len(C))

# Requests: take all r values used in pickup/dropoff  
r_values = nodes[np.isin(nodes[:,1], ["pickup","dropoff","transfer"]), 2]
r_values = r_values[r_values!=""]
R = sorted(set(r_values.astype(int)))

# Vehicles
V = 2   # fleet size
K = list(range(V))  # vehicle set with (future) travel nodes for each vehicle

# Depots
zeroDepot = int(nodes[nodes[:,1]=="depot0",0][0])
endDepot  = int(nodes[nodes[:,1]=="depot1",0][0])

# C_r sets (transfer nodes per request)
Cr = {r: nodes[(nodes[:,1]=="transfer") & (nodes[:,2]==str(r)),0].astype(int).tolist() for r in R}

##### Parametres #####

# Service durations, time windows, load
di = {int(row[0]): float(row[5]) for row in nodes if row[5] != ""}
ei = {int(row[0]): float(row[3]) for row in nodes if row[3] != ""}
li = {int(row[0]): float(row[4]) for row in nodes if row[4] != ""}
qr = {int(r): 1 for r in R}   # simple: each request has load 1 
## qr = {int(row[0]): float(row[6]) for row in nodes if row[6] != ""}

# Travel time and cost dictionaries
tij = {(i,j): float(t_transfer[i,j]) for i in range(len(t_transfer)) for j in range(len(t_transfer))}
# print(tij[0,11])

# Example planning horizon etc.
Q = 4
T = 3600*4
ek = {k: 0 for k in range(len(K))}
lk = {k: 86400 for k in range(len(K))}
# fi_r = {(r,i): 1 if i in P and nodes[nodes[:,0]==i,2][0]==str(r) else -1 if i in D and nodes[nodes[:,0]==i,2][0]==str(r) else 0 for r in R for i in N}
fi_r = {}
for r in R:
    for i in N:  
        fi_r[r,i] = 0
     
for r in R:
    fi_r[r,r] = 1
    fi_r[r,r+len(R)] = -1

Lbar = {i: 1800 for i in P}
pair_pi_di = {p:d for p,d in zip(P,D)}
M = 100000
n=len(P)  # number of requests
n_K=len(K)  # number of vehicles

# ei and li bounds for dropoff nodes
for key in pair_pi_di.keys():
    i = pair_pi_di[key]
    # print(i)
    e = ei[key] + di[key] + tij[key, i]
    l = li[key] + di[key] + tij[key, i]
    ei[i] = e
    li[i] = l
# print(ei)
# print(li)

# Electric Vehicle Parametres
C_bat_kWh = 100 # Battery Capacity in kWh
C_bat = C_bat_kWh*(60*60) # Battery Capacity in kWs
alpha = 200 # Charging Power in kW
beta = 15 # Average Power in kW
gamma = 0.1 # Minimum battery capacity set at 10%

tij_original = tij.copy()

##### Arc Elimination #####
for i in N:   # remove arcs to/from depot for non-depot nodes
    if (i,0) in tij:
        del tij[i,0]   # Remove arcs to depot0 from any node except depot0
    if (2*n+1,i) in tij:
        del tij[2*n+1,i]  # Remove arcs from depot1 to any node except depot1

for i in D:
    if (0,i) in tij:
        del tij[0,i]   # Remove arcs to dropoff from depot0

for i in P:
    if (i,2*n+1) in tij:
        del tij[i,2*n+1]   # Remove arcs to depot1 from any pickup

for i in N:
    if (i,i) in tij:
        del tij[i,i]   # Remove loops

for i in P:
    if (n+i,i) in tij:
        del tij[n+i,i]   # Remove arcs from drop-off to its pickup

# Remove infeasible arcs due to time windows
for i in N:
    for j in N:
        if (i, j) in tij:
            # Check if ei[i], di[i], and li[j] exist
            if i in ei and i in di and j in li:
                if ei[i] + di[i] + tij[i, j] >= li[j]:
                    del tij[i, j]

for i in P:
    for j in N:
        if (i,j) in tij and (j,n+i) in tij:
            if tij[i,j] + di[j] + tij[j,n+i] >= Lbar[i]:
                del tij[i,j]
                del tij[j,n+i]
                # Remove infeasible arcs due to ride time

# Path infeasibility: if {j, i, n+j, n+i} infeasible, remove (i, n+j)
for i in P:
    for j in P:
        if i != j:
            if (j, i) in tij and (i, n+j) in tij and (n+j, n+i) in tij:
                if tij[j, i] + di[i] + tij[i, n+j] + di[n+j] + tij[n+j, n+i] > Lbar[i]:
                    del tij[i, n+j]

# Path infeasibility: if {i, n+i, j, n+j} infeasible, remove (n+i, j)
for i in P:
    for j in P:
        if i != j:
            if (i, n+i) in tij and (n+i, j) in tij and (j, n+j) in tij:
                if tij[i, n+i] + di[n+i] + tij[n+i, j] + di[j] + tij[j, n+j] > Lbar[i]:
                    del tij[n+i, j]

for Si in C:
    if (Si,Si) in tij:
        del tij[Si,Si]   # Remove loops at transfer nodes

# Correct transfer node handling
for i in P:
    drop_node = n + i
    for Si in Cr[i]:
        if (drop_node, Si) in tij:
            del tij[drop_node, Si]  # Remove arcs from drop-off to transfer node for same request
        if (Si, i) in tij:
            del tij[Si, i]  # Remove arcs from transfer to pick-up node for same request

for i in P:
    for j in P:
        if i != j and (i, n+j) not in tij:
            for Si in Cr[i]:
                for Sj in Cr[j]:
                    if (Si, Sj) in tij:
                        del tij[Si, Sj]
                if (Si, n+j) in tij:
                    del tij[Si, n+j]
                if (Si, j) in tij:
                    del tij[Si, j]

cij = tij.copy()

##### Export Data #####
data = dict(nodes=nodes, N=N, P=P, D=D, C=C, R=R, K=K, Cr=Cr, F=F,
            Q=Q, T=T, cij=cij, tij=tij, tij_original=tij_original, qr=qr, di=di, ei=ei, li=li,
            ek=ek, lk=lk, fi_r=fi_r, Lbar=Lbar, endDepot=endDepot,
            zeroDepot=zeroDepot, pair_pi_di=pair_pi_di, M=M, n=n, n_K=n_K, C_bat=C_bat, alpha=alpha, beta=beta, gamma=gamma)

