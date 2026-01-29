import argparse
from pathlib import Path
import pyodbc
from pathlib import Path

BASE_DIR = Path(r"C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes"
                r"\\Stage-Inria-01-09-2025--28-02-2026\\code_and_instances")

SERVER = r"localhost\SQLEXPRESS"
DATABASE = "darp"

data_darp = ["datafile0.txt", "datafile1.txt", "datafile2.txt", "datafile3.txt", "datafile4.txt", "datafile5.txt", "datafile6.txt", "datafile7.txt", "datafile8.txt", "datafile9.txt" ]

SCHEMA_SQL = """
IF OBJECT_ID('dbo.instances', 'U') IS NULL
BEGIN
    CREATE TABLE dbo.instances (
        id INT IDENTITY(1,1) PRIMARY KEY,
        name NVARCHAR(255) NOT NULL UNIQUE,
        source_path NVARCHAR(1024) NOT NULL,
        n_vehicles INT NOT NULL,
        n_requests INT NOT NULL,
        n_terminals INT NOT NULL,
        created_at DATETIME2 NOT NULL DEFAULT SYSUTCDATETIME()
    );
END;

IF OBJECT_ID('dbo.nodes', 'U') IS NULL
BEGIN
    CREATE TABLE dbo.nodes (
        id INT IDENTITY(1,1) PRIMARY KEY,
        instance_id INT NOT NULL,
        node_id INT NOT NULL,
        x FLOAT NOT NULL,
        y FLOAT NOT NULL,
        twL FLOAT NOT NULL,
        twU FLOAT NOT NULL,
        mobile INT NOT NULL,
        wheelchair INT NOT NULL,
        service_dur FLOAT NOT NULL,
        CONSTRAINT FK_nodes_instances FOREIGN KEY (instance_id) REFERENCES dbo.instances(id) ON DELETE CASCADE,
        CONSTRAINT UQ_nodes_instance_node UNIQUE (instance_id, node_id)
    );

    CREATE INDEX IX_nodes_instance_id ON dbo.nodes(instance_id);
END;
"""


def parse_header(line: str) -> tuple[int, int, int]:
    parts = line.strip().split()
    if len(parts) < 3:
        raise ValueError(f"Bad header line (need 3 ints): {line!r}")
    return int(parts[0]), int(parts[1]), int(parts[2])


def parse_row(line: str) -> tuple[int, float, float, float, float, int, int, float]:
    parts = line.strip().split()
    if len(parts) != 8:
        raise ValueError(f"Bad row (need 8 fields): {line!r}")
    return (
        int(parts[0]),
        float(parts[1]),
        float(parts[2]),
        float(parts[3]),
        float(parts[4]),
        int(parts[5]),
        int(parts[6]),
        float(parts[7]),
    )


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


def import_txt_to_sql_express(txt_path: Path, conn: pyodbc.Connection, file: str, replace: bool = False) -> int:
    if not txt_path.exists():
        raise FileNotFoundError(txt_path)

    cur = conn.cursor()
    # cur.execute(SCHEMA_SQL)

    instance_name = txt_path.stem

    if replace:
        cur.execute("SELECT id FROM dbo.instances WHERE name = ?", instance_name)
        row = cur.fetchone()
        if row:
            cur.execute("DELETE FROM dbo.instances WHERE id = ?", int(row[0]))

    with txt_path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        if not header:
            raise ValueError("File is empty")

        n_vehicles, n_requests, n_terminals = parse_header(header)

        cur.execute(
            """
            INSERT INTO dbo.instances (name, source_path, n_vehicles, n_requests, n_terminals)
            OUTPUT INSERTED.id
            VALUES (?, ?, ?, ?, ?)
            """,
            instance_name, str(file), n_vehicles, n_requests, n_terminals
        )
        instance_id = int(cur.fetchone()[0])

        batch = []
        line_no = 1
        for line in f:
            line_no += 1
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            try:
                batch.append(parse_row(s))
            except Exception as e:
                raise ValueError(f"Parse error line {line_no}: {e}") from e

        # Fast insert
        cur.fast_executemany = True
        cur.executemany(
            """
            INSERT INTO dbo.nodes (instance_id, node_id, x, y, twL, twU, mobile, wheelchair, service_dur)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [(instance_id, *r) for r in batch]
        )

    conn.commit()
    return instance_id


if __name__ == "__main__":

    with connect_sql_express(server=SERVER, database=DATABASE, trusted=True) as conn:
        for file in data_darp:
            txt_path = BASE_DIR / file
            iid = import_txt_to_sql_express(txt_path, conn, file, replace=True)
            print(f"OK. Imported {file} → instance_id={iid}")
