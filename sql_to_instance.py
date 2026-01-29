import argparse
from pathlib import Path
import pyodbc
import numpy as np

BASE_DIR = Path(r"C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes"
                r"\\Stage-Inria-01-09-2025--28-02-2026\\code_and_instances")

SERVER = r"localhost\SQLEXPRESS"
DATABASE = "darp"

def connect_sql_express(server: str, database: str) -> pyodbc.Connection:

    conn_str = (
        "Driver={ODBC Driver 18 for SQL Server};"
        f"Server={server};"
        f"Database={database};"
        "Trusted_Connection=yes;"
        "Encrypt=yes;"
        "TrustServerCertificate=yes;"
    )

    return pyodbc.connect(conn_str, autocommit=False)


def select_earliest_n_and_offset_nodes(conn: pyodbc.Connection, n, offset, instance_id):

    sql = """SELECT node_id, instance_id, x, y, node_type, twL_constrained, twU_constrained, mobile, wheelchair, service_dur
            FROM dbo.nodes
            WHERE instance_id = ? AND node_type = 'Pickup Node'
            ORDER BY twL_constrained
            OFFSET ? ROWS
            FETCH NEXT ? ROWS ONLY;"""
    
    cur = conn.cursor()

    cur.execute(sql, instance_id, offset, n)

    new_request_nodes = cur.fetchall()

    return new_request_nodes


def select_n_random_requests(conn: pyodbc.Connection, n, instance_id):

    sql_pickup = """SELECT TOP (?) node_id, instance_id, x, y, node_type,
        twL, twU,
        mobile, wheelchair, service_dur
        FROM dbo.nodes
        WHERE instance_id = ? AND node_type = 'Pickup Node'
        ORDER BY NEWID();"""

    cur1 = conn.cursor()
    cur2 = conn.cursor()

    cur1.execute(sql_pickup, n, instance_id)

    new_pickup_nodes = cur1.fetchall()

    pickup_node_ids = [new_pickup_nodes[i][0] for i in range(n)]

    sql_dropoff = """SELECT node_id, instance_id, x, y, node_type,
        twL, twU,
        mobile, wheelchair, service_dur
        FROM dbo.nodes
        WHERE instance_id = ? AND node_id = ?"""
    
    new_dropoff_nodes = []

    for pickup_node_id in pickup_node_ids:
    
        cur2.execute(sql_dropoff, instance_id, pickup_node_id + 200)

        new_dropoff_node = cur2.fetchone()

        new_dropoff_nodes.append(new_dropoff_node)

    return new_pickup_nodes, new_dropoff_nodes


def select_transfer_nodes(conn, instance_id):

    sql = """SELECT node_id, instance_id, x, y, node_type, twL_constrained, twU_constrained, mobile, wheelchair, service_dur
        FROM dbo.nodes
        WHERE instance_id = ? 
            AND node_type = 'Transfer Node'
            AND x IN (-20, -10, 0, 10, 20)
            AND y IN (-20, -10, 0, 10, 20);"""
    
    cur = conn.cursor()

    cur.execute(sql, instance_id)

    new_request_nodes = cur.fetchall()

    return new_request_nodes





def build_pickup_node(idx, pickup_node):
    ## np.array(Node id, Node type, request, Lower Time Window, Upper Time Window, Service Duration, Passenger Capacity)
    node_id = idx
    request = idx
    node_type = "pickup"
    twL = pickup_node[5]
    twU = pickup_node[6]
    di = pickup_node[9]
    qr = pickup_node[7] + pickup_node[8]
    x = pickup_node[2]
    y = pickup_node[3]

    node_format = [node_id, node_type, request, twL, twU, di, qr, float(x), float(y)]

    return node_format

def build_dropoff_node(idx, dropoff_node):
    node_id = idx 
    request = idx - 10
    node_type = "dropoff"
    twL = dropoff_node[5]
    twU = dropoff_node[6]
    di = dropoff_node[9]
    qr = 0
    x = dropoff_node[2]
    y = dropoff_node[3]

    node_format = [node_id, node_type, request, twL, twU, di, qr, float(x), float(y)]
    
    return node_format

def build_transfer_node(idx, transfer_node):
    node_id = idx 
    request = ''
    node_type = "transfer"
    twL = 0
    twU = 1440
    di = transfer_node[9]
    qr = transfer_node[7] + transfer_node[8]
    x = transfer_node[2]
    y = transfer_node[3]

    node_format = [node_id, node_type, request, twL, twU, di, qr, x, y]
    
    return node_format

def build_nodes(new_pickup_nodes, new_dropoff_nodes, new_transfer_nodes):
    nodes=[]
    depot_node = [0, "depot0", "", 0, 1440, 0, 0]
    idx = 0
    nodes.append([idx, depot_node[1], depot_node[2], depot_node[3], depot_node[4], depot_node[5], depot_node[6], 0, 0])
    for pickup_node in new_pickup_nodes:
        idx += 1
        nodes.append(build_pickup_node(idx, pickup_node))
    
    for dropoff_node in new_dropoff_nodes:
        idx += 1
        nodes.append(build_dropoff_node(idx, dropoff_node))

    idx += 1
    nodes.append([idx, "depot1", depot_node[2], depot_node[3], depot_node[4], depot_node[5], depot_node[6], 0, 0])

    for transfer_node in new_transfer_nodes:
        idx += 1
        nodes.append(build_transfer_node(idx, transfer_node))

    return nodes


def dist_i_j(node_i,node_j):
    xi = node_i[7]
    xj = node_j[7]
    yi = node_i[8]
    yj = node_j[8]
    return np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)


def pt_terminal_distance(node_i, node_j, speed_factor_pt: float = 1.0, big_m: float = 10000.0) -> float:
    xi = node_i[7]
    xj = node_j[7]
    yi = node_i[8]
    yj = node_j[8]

    if node_i == node_j:
        return 0

    # Exclude trips involving more than one change (Molenbruch et al. logic)
    if ((xi == -xj and abs(yi) != 10 and abs(yj) != 10) or
        (yi == -yj and abs(xi) != 10 and abs(xj) != 10)):
        return big_m

    # Manhattan travel time, then scaled
    manhattan = abs(xi - xj) + abs(yi - yj)
    return manhattan * speed_factor_pt


def n_closest_transfer(n, node, transfer_nodes):
    closest = [(tn[0], dist_i_j(node, tn)) for tn in transfer_nodes]
    closest.sort(key=lambda x: x[1])
    return closest[:n]


def automatic_dist_i_j(node_i, node_j, transfer_nodes, n_closest=3, speed_factor_pt=0.5, big_m=10000.0):
    i_type = node_i[1]
    j_type = node_j[1]

    # transfer <-> transfer uses PT metric
    if i_type == 'transfer' and j_type == 'transfer':
        return pt_terminal_distance(node_i, node_j, speed_factor_pt, big_m)

    # non-transfer <-> transfer: only allow if transfer is among n closest to the non-transfer node
    if i_type == 'transfer' and j_type != 'transfer':
        allowed_ids = {tid for tid, _ in n_closest_transfer(n_closest, node_j, transfer_nodes)}
        return dist_i_j(node_i, node_j) if node_i[0] in allowed_ids else big_m

    if i_type != 'transfer' and j_type == 'transfer':
        allowed_ids = {tid for tid, _ in n_closest_transfer(n_closest, node_i, transfer_nodes)}
        return dist_i_j(node_i, node_j) if node_j[0] in allowed_ids else big_m

    # otherwise: standard euclidean
    return dist_i_j(node_i, node_j)


def build_tij(nodes):
    len_nodes = len(nodes)
    tij = np.zeros((len_nodes, len_nodes))
    transfer_nodes = [node for node in nodes if node[1] == 'transfer']
    for i in range(len_nodes):
        for j in range(len_nodes):
            tij[i][j] = round(automatic_dist_i_j(nodes[i],nodes[j], transfer_nodes), 1)
    return tij


def build_random_instance():
    n = 5
    instance_id = 1
    conn = connect_sql_express(".\\SQLEXPRESS", "darp")
    new_pickup_nodes, new_dropoff_nodes = select_n_random_requests(conn, n, instance_id = instance_id)
    new_transfer_nodes = select_transfer_nodes(conn, instance_id = instance_id)
    nodes = build_nodes(new_pickup_nodes, new_dropoff_nodes, new_transfer_nodes)
    tij = build_tij(nodes)
    nodes = [row[:-2] for row in nodes]
    return nodes, tij



if __name__ == "__main__":
    nodes, tij = build_random_instance()
    print(nodes)
    # print(new_pickup_nodes)
    # print("new_dropoff_nodes is: \n")
    # print(new_dropoff_nodes)
    # print("new_transfer_nodes is: \n")
    # print(new_transfer_nodes)


    # print("New Nodes are: \n")
    # for node in nodes:
    #     print(node)


    # nodes = nodes[:, :-2]
    # print(nodes)
    # print(tij)
    # print(tij[1])
    # print("new_request_nodes is: \n")
    # print(new_dropoff_nodes)
    # for i in range(n):
    #     print(f"Request {i} is :")
    #     print(new_pickup_nodes[i])
    #     print(new_dropoff_nodes[i], "\n")
    # for transfer_node in new_transfer_nodes:
    #     print("Transfer nodes are: ", transfer_node)



















