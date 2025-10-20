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
    [10,  "transfer","", "" ,   "",    0,   0],
    [11,  "transfer","", "" ,   "",    0,   0],
    [12,  "transfer","", "" ,   "",    0,   0],
], dtype=object)

# Example: build 3 physical transfer stops × 4 requests -> ids 10..21
# print("Nodes:\n", nodes)

##### Build sets #####
max_visits_transfer = 10
max_visits_node = 2

# Nodes
N_minor = nodes[:,0].astype(int).tolist()
# N_type = nodes[:,1].astype(str).tolist()
P_minor = nodes[nodes[:,1]=="pickup",0].astype(int).tolist()
D_minor = nodes[nodes[:,1]=="dropoff",0].astype(int).tolist()
C_minor = nodes[nodes[:,1]=="transfer",0].astype(int).tolist()
# Depots
zeroDepot = int(nodes[nodes[:,1]=="depot0",0][0])
endDepot  = int(nodes[nodes[:,1]=="depot1",0][0])
zeroDepot_node = (zeroDepot,1)
endDepot_node = (endDepot,1)

P = [(i,1) for i in P_minor]
D = [(i,1) for i in D_minor]
C = [(i,j) for j in range(max_visits_transfer) for i in C_minor]
N = [zeroDepot_node, endDepot_node] + P + D + C

# print("P is: ",P)
# print("D is: ",D)
# print("C is: ",C)
# print("N is: ",N)

# print(len(N), len(P), len(D), len(C))

# Departure every 5 minutes and 30 seconds stop time
# D_10_11 = [300 * i for i in range(24*60*60/300)]
# D_11_12 = [300 * i + t[10,11] + 30 for i in range(24*60*60/300)]
# D_10_12 = [300 * i for i in range(24*60*60/300)]
# D_11_10 = [300 * i + t[12,11] + 30 for i in range(24*60*60/300)]

Departures = {}

interval = 300            
planning_horizon = 24 * 60 * 60
n_intervals = planning_horizon // interval

# Transfer stops in order
transfer_nodes = sorted(C_minor)

# === Forward direction (left → right) ===
left_terminal = transfer_nodes[0]
for j in range(1, len(transfer_nodes)):
    Departures[(left_terminal, transfer_nodes[j])] = {a:interval * a for a in range(int(n_intervals))}

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
    Departures[(right_terminal, transfer_nodes[j])] ={a:interval * a for a in range(int(n_intervals))}

# Backward arcs (including skips)
for i in range(len(transfer_nodes)-2, -1, -1):
    for j in range(i-1, -1, -1):
        Departures[(transfer_nodes[i], transfer_nodes[j])] = {
            a: float(d + sum(t[transfer_nodes[k], transfer_nodes[k-1]] for k in range(i, j, -1)))
            for a, d in Departures[(right_terminal, transfer_nodes[i])].items()
        }

# print("Departures:", Departures)


           
# Requests: take all r values used in pickup/dropoff  
r_values = nodes[np.isin(nodes[:,1], ["pickup","dropoff"]), 2]
r_values = r_values[r_values != ""]
R = sorted(set(r_values.astype(int)))

# Vehicles
V = 2   # fleet size
K = list(range(V))  # vehicle set with (future) travel nodes for each vehicle

##### Parametres #####
  
# Service durations, time windows, load
di = {int(row[0]): float(row[5]) for row in nodes if row[5] != ""}
ei = {int(row[0]): float(row[3]) for row in nodes if row[3] != ""}
li = {int(row[0]): float(row[4]) for row in nodes if row[4] != ""}
qr = {int(row[0]): float(row[6]) for row in nodes if row[6] != ""}
# print(di)
# print(ei)
# print(li)
# print(qr)

# Travel time and cost dictionaries
tij = {(i,j): float(t[i,j]) for i in range(len(t)) for j in range(len(t))}

# Example planning horizon etc.
Q = 4
T = 3600*4
ek = {k: 0 for k in range(len(K))}
lk = {k: 86400 for k in range(len(K))}
fi_r = {}
for r in R:
    for i in N:
        fi_r[r,i] = 0

for r in R:
    fi_r[r,r] = 1
    fi_r[r,r+len(R)] = -1

Lbar = {i[0]: 1800 for i in P}
# for p,d in zip(P,D):
#     print("p,d is :",p,d)
pair_pi_di = {p:d for p,d in zip(P,D)}
# print("pair_pi_di is :", pair_pi_di)
M = 100000
n=len(P)  # number of requests
n_K=len(K)  # number of vehicles

# ei and li bounds for dropoff nodes
for key in pair_pi_di.keys():
    i = pair_pi_di[key][0]
    # print(key)
    # print(pair_pi_di[key])
    # print(i)
    e = ei[key[0]] + di[key[0]] + tij[key[0], pair_pi_di[key][0]]
    l = li[key[0]] + di[key[0]] + tij[key[0], pair_pi_di[key][0]]
    ei[i] = e
    li[i] = l
# print(ei)
# print(li)

##### Arc Elimination #####
for i in N_minor:   # remove arcs to/from depot for non-depot nodes
    if (i,zeroDepot) in tij:
        del tij[i,zeroDepot]   # Remove arcs to depot0 from any node except depot0
    if (endDepot,i) in tij:
        del tij[endDepot,i]  # Remove arcs from depot1 to any node except depot1

for i in D_minor:
    if (zeroDepot,i) in tij:
        del tij[zeroDepot,i]   # Remove arcs to dropoff from depot0

for i in P_minor:
    if (i,endDepot) in tij:
        del tij[i,endDepot]   # Remove arcs to depot1 from any pickup

for i in N_minor:
    if (i,i) in tij:
        del tij[i,i]   # Remove loops

for i in P_minor:
    if (n+i,i) in tij:
        del tij[n+i,i]   # Remove arcs from drop-off to its pickup

# Remove infeasible arcs due to time windows
for i in N_minor:
    for j in N_minor:
        if (i, j) in tij:
            # Check if ei[i], di[i], and li[j] exist
            if i in ei and i in di and j in li:
                if ei[i] + di[i] + tij[i, j] >= li[j]:
                    del tij[i, j]

for i in P_minor:
    for j in N_minor:
        if (i,j) in tij and (j,n+i) in tij:
            if tij[i,j] + di[j] + tij[j,n+i] >= Lbar[i]:
                del tij[i,j]
                del tij[j,n+i]
                # Remove infeasible arcs due to ride time

# Path infeasibility: if {j, i, n+j, n+i} infeasible, remove (i, n+j)
for i in P_minor:
    for j in P_minor:
        if i != j:
            if (j, i) in tij and (i, n+j) in tij and (n+j, n+i) in tij:
                if tij[j, i] + di[i] + tij[i, n+j] + di[n+j] + tij[n+j, n+i] > Lbar[i]:
                    del tij[i, n+j]

# Path infeasibility: if {i, n+i, j, n+j} infeasible, remove (n+i, j)
for i in P_minor:
    for j in P_minor:
        if i != j:
            if (i, n+i) in tij and (n+i, j) in tij and (j, n+j) in tij:
                if tij[i, n+i] + di[n+i] + tij[n+i, j] + di[j] + tij[j, n+j] > Lbar[i]:
                    del tij[n+i, j]

for Si in C_minor:
    if (Si,Si) in tij:
        del tij[Si,Si]   # Remove loops at transfer nodes

# print(nodes[nodes[:, 0] == 10, 1][0])

cij = tij.copy()

##### Export Data #####
data = dict(nodes=nodes, N=N, P=P, D=D, C=C, R=R, K=K, N_minor=N_minor, P_minor=P_minor, D_minor=D_minor, C_minor=C_minor, zeroDepot_node=zeroDepot_node, endDepot_node=endDepot_node, Departures=Departures,
            Q=Q, T=T, cij=cij, tij=tij, qr=qr, di=di, ei=ei, li=li,
            ek=ek, lk=lk, fi_r=fi_r, Lbar=Lbar, endDepot=endDepot,
            zeroDepot=zeroDepot, pair_pi_di=pair_pi_di, M=M, n=n, n_K=n_K)
