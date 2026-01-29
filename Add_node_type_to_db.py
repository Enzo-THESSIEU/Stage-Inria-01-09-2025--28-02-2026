import pyodbc

def produce_conn_str(SERVER, DATABASE):
    conn_str = (
            "Driver={ODBC Driver 18 for SQL Server};"
            f"Server={SERVER};"
            f"Database={DATABASE};"
            "Trusted_Connection=yes;"
            "Encrypt=no;"
        )
    return conn_str


def add_column_to_table(conn_str, TABLE, NEW_COLUMN, COLUMN_TYPE):
    sql = f"""
        IF COL_LENGTH('{TABLE}', '{NEW_COLUMN}') IS NULL
        BEGIN
            ALTER TABLE {TABLE}
            ADD {NEW_COLUMN} {COLUMN_TYPE} NULL;
        END
        """
    
    with pyodbc.connect(conn_str) as conn:
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()

    print(f"Column '{NEW_COLUMN}' added to {TABLE} (if it did not already exist)")


def extract_instance_ids(conn_str):
    TABLE = "dbo.instances"
    sql = f"""
    SELECT id
    FROM {TABLE}
    ORDER BY id;
    """
    with pyodbc.connect(conn_str) as conn:
        cur = conn.cursor()
        cur.execute(sql)
        instances = cur.fetchall()

    print("instances is: ", instances)

    return [instance.id for instance in instances]

    
def extract_instance_parametres(INSTANCE_ID, conn_str):
    TABLE = "dbo.instances"
    sql = f"""
    SELECT n_vehicles, n_requests, n_terminals
    FROM {TABLE}
    WHERE id = ?;
    """

    with pyodbc.connect(conn_str) as conn:
        cur = conn.cursor()
        cur.execute(sql, (INSTANCE_ID))
        n_values = cur.fetchall()

    print("n_values is: ", n_values)

    return n_values[0][0], n_values[0][1], n_values[0][2]


def define_node_type(INSTANCE_ID, conn_str):
    n_vehicles, n_requests, n_terminals = extract_instance_parametres(INSTANCE_ID, conn_str)
    pickup_min = 0
    pickup_max = n_requests - 1
    dropoff_min = pickup_max + 1
    dropoff_max = pickup_max + n_requests 
    depot_min = dropoff_max + 1
    depot_max = dropoff_max + n_vehicles 
    transfer_min = depot_max + 1
    transfer_max = depot_max + n_terminals 

 
    sql = """
        UPDATE dbo.nodes
        SET node_type =
            CASE
                WHEN node_id BETWEEN ? AND ? THEN 'Pickup Node'
                WHEN node_id BETWEEN ? AND ? THEN 'Dropoff Node'
                WHEN node_id BETWEEN ? AND ? THEN 'Depot Node'
                WHEN node_id BETWEEN ? AND ? THEN 'Transfer Node'
                ELSE 'Other'
            END
        WHERE instance_id = ?;
        """

    with pyodbc.connect(conn_str) as conn:
        cur = conn.cursor()
        cur.execute(sql, (
            pickup_min, pickup_max,
            dropoff_min, dropoff_max,
            depot_min, depot_max,
            transfer_min, transfer_max,
            INSTANCE_ID
        ))
        conn.commit()

    print("Node Type added for instance: ", INSTANCE_ID)


if __name__ == "__main__":
    # ---- connection ----
    SERVER = r".\SQLEXPRESS"          # or r"(localdb)\MSSQLLocalDB"
    DATABASE = "darp"                # choose DB here

    conn_str = produce_conn_str(SERVER, DATABASE)

    add_column_to_table(conn_str, "dbo.nodes", "node_type", "VARCHAR(50)")

    instances = extract_instance_ids(conn_str)
    for INSTANCE_ID in instances:
        define_node_type(INSTANCE_ID, conn_str)

    print("All node tpes added to SQL Table for each instance")


