import ast
import pandas as pd
from pyvis.network import Network
import math

# --- Node mapping for duplicates ---
NODE_MAP = {
    13: 10, 16: 10, 19: 10,   # duplicates of 10
    14: 11, 17: 11, 20: 11,   # duplicates of 11
    15: 12, 18: 12, 21: 12    # duplicates of 12
}

def remap_node(node):
    """Return remapped node if it is a duplicate; otherwise return itself."""
    return NODE_MAP.get(node, node)


def parse_vehicle_routes(row):
    """Parse all vehicle route columns from the CSV row and return a dict."""
    vehicle_routes = {}
    for col in row.index:
        if "Vehicle" in col and "Route" in col and isinstance(row[col], str) and row[col].strip():
            try:
                data = ast.literal_eval(row[col])
                arcs = []
                for elem in data:
                    if isinstance(elem, list) and len(elem) >= 1:
                        # First element is typically [from, to]
                        if isinstance(elem[0], list) and len(elem[0]) == 2:
                            i, j = elem[0]
                        elif len(elem) >= 2 and isinstance(elem[0], int) and isinstance(elem[1], int):
                            i, j = elem[0], elem[1]
                        else:
                            continue

                        # --- Apply remapping ---
                        i_mapped = remap_node(i)
                        j_mapped = remap_node(j)
                        arcs.append((i_mapped, j_mapped))
                        # print(f"Arc ({i},{j}) → ({i_mapped},{j_mapped})")

                veh_num = int(col.split()[1])
                vehicle_routes[veh_num] = arcs
            except Exception as e:
                print(f"⚠️ Error parsing {col}: {e}")
        # print("vehicle_routes:", vehicle_routes)
    return vehicle_routes


def display_vehicle_routes_from_idarp(csv_path, model_name=None, title="IDARP Vehicle Routes Visualization"):
    """
    Visualize vehicle routes from an IDARP model CSV file.
    Automatically detects columns like 'Vehicle 1 Route', 'Vehicle 2 Route', etc.
    """

    # --- Load CSV ---
    df = pd.read_csv(csv_path)

    # If user didn't specify which model, take the first one
    if model_name:
        row = df[df["Model"] == model_name].iloc[0]
    else:
        row = df.iloc[0]
        model_name = row["Model"]

    print(f"🔍 Visualizing model: {model_name}")

    # --- Parse vehicle routes ---
    vehicle_routes = parse_vehicle_routes(row)

    # --- Create nodes (0–12 on a circle) ---
    nodes = {}
    for node_id in range(13):
        angle = 2 * math.pi * node_id / 13
        x = math.cos(angle)
        y = math.sin(angle)
        nodes[node_id] = (x, y, f"Node {node_id}")

    # --- Collect all possible arcs for faint background ---
    arcs = [(i, j) for i in nodes for j in nodes if i != j]

    # --- Create PyVis network ---
    net = Network(
        height='750px',
        width='100%',
        bgcolor='#20232a',
        font_color='white',
        directed=True
    )
    net.barnes_hut(gravity=-80000, central_gravity=0.3)

# Add nodes with grouped colors
    for node_id, (x, y, label) in nodes.items():
        if node_id in [0, 9]:
            color = "#4CAF50"  # Green for depots or special nodes
        elif node_id in [1, 2, 3, 4]:
            color = "#e67e22"  # Orange for first group
        elif node_id in [5, 6, 7, 8]:
            color = "#9b59b6"  # Purple for second group
        else:
            color = "#3498db"  # Blue for Transfer Nodes

        net.add_node(
            node_id,
            label=label,
            x=x * 300,
            y=y * 300,
            physics=False,
            color=color
        )



    # Faint gray arcs for full connectivity
    for i, j in arcs:
        net.add_edge(i, j, color="rgba(200,200,200,0.05)", width=1)

    # --- Draw vehicle arcs ---
    colors = ["#e74c3c", "#3498db", "#f1c40f", "#9b59b6", "#1abc9c", "#e67e22", "#2ecc71"]
    for idx, (veh, route) in enumerate(vehicle_routes.items()):
        color = colors[idx % len(colors)]
        for (i, j) in route:
            if i in nodes and j in nodes:
                net.add_edge(i, j, color=color, width=4, title=f"Vehicle {veh}")

    # --- Save output ---
    output_file = f"{model_name.replace(' ', '_').lower()}_routes.html"
    net.write_html(output_file)
    print(f"✅ Graph written to {output_file}")


# Example usage
if __name__ == "__main__":
    csv_path = r"C:\Users\enzot\Documents\Césure\1ère césure inria Lille\Codes\results_full.csv"
    display_vehicle_routes_from_idarp(csv_path, model_name="IDARP_Model_dup1_arc1_var1_sub1_tim0_use0")


# import ast
# import pandas as pd
# from pyvis.network import Network
# import math

# def parse_vehicle_routes(row):
#     """Parse all vehicle route columns from the CSV row and return a dict."""
#     vehicle_routes = {}
#     for col in row.index:
#         if "Vehicle" in col and "Route" in col and isinstance(row[col], str) and row[col].strip():
#             try:
#                 data = ast.literal_eval(row[col])
#                 arcs = []
#                 for elem in data:
#                     if isinstance(elem, list) and len(elem) >= 1:
#                         # First element is typically [from, to]
#                         if isinstance(elem[0], list) and len(elem[0]) == 2:
#                             arcs.append(tuple(elem[0]))
#                         elif len(elem) >= 2 and isinstance(elem[0], int) and isinstance(elem[1], int):
#                             arcs.append((elem[0], elem[1]))
#                 veh_num = int(col.split()[1])
#                 vehicle_routes[veh_num] = arcs
#             except Exception as e:
#                 print(f"⚠️ Error parsing {col}: {e}")
#     return vehicle_routes


# def parse_request_and_transfer_paths(row):
#     """Try to extract request and transfer paths if available in the CSV."""
#     request_paths = []
#     transfer_paths = []

#     for col in row.index:
#         val = row[col]
#         if not isinstance(val, str) or not val.strip():
#             continue
#         try:
#             data = ast.literal_eval(val)
#             if "request" in col.lower():
#                 # Expect list of (pickup, dropoff)
#                 for elem in data:
#                     if isinstance(elem, list) and len(elem) == 2:
#                         request_paths.append(tuple(elem))
#             elif "transfer" in col.lower():
#                 # Expect list of arcs (i,j) for transfer nodes
#                 for elem in data:
#                     if isinstance(elem, list) and len(elem) == 2:
#                         transfer_paths.append(tuple(elem))
#         except Exception:
#             continue

#     return request_paths, transfer_paths


# def display_vehicle_routes_from_idarp(csv_path, model_name=None, title="IDARP Vehicle Routes Visualization"):
#     """
#     Visualize vehicle routes from an IDARP model CSV file.
#     Adds request and transfer paths if available.
#     """

#     df = pd.read_csv(csv_path)

#     if model_name:
#         row = df[df["Model"] == model_name].iloc[0]
#     else:
#         row = df.iloc[0]
#         model_name = row["Model"]

#     print(f"🔍 Visualizing model: {model_name}")

#     # --- Parse vehicle, request, and transfer data ---
#     vehicle_routes = parse_vehicle_routes(row)
#     request_paths, transfer_paths = parse_request_and_transfer_paths(row)

#     # --- Create node positions (0–12 on a circle) ---
#     nodes = {}
#     for node_id in range(13):
#         angle = 2 * math.pi * node_id / 13
#         x = math.cos(angle)
#         y = math.sin(angle)
#         nodes[node_id] = (x, y, f"Node {node_id}")

#     arcs = [(i, j) for i in nodes for j in nodes if i != j]

#     # --- Create PyVis network ---
#     net = Network(
#         height='750px',
#         width='100%',
#         bgcolor='#20232a',
#         font_color='white',
#         directed=True
#     )
#     net.barnes_hut(gravity=-80000, central_gravity=0.3)

#     # Add nodes
#     for node_id, (x, y, label) in nodes.items():
#         net.add_node(
#             node_id,
#             label=label,
#             x=x * 300,
#             y=y * 300,
#             physics=False,
#             color="#4CAF50" if node_id == 0 else "#3498db"
#         )

#     # Background faint arcs
#     for i, j in arcs:
#         net.add_edge(i, j, color="rgba(200,200,200,0.05)", width=1)

#     # --- Draw vehicle arcs ---
#     colors = ["#e74c3c", "#3498db", "#f1c40f", "#9b59b6", "#1abc9c", "#e67e22", "#2ecc71"]
#     for idx, (veh, route) in enumerate(vehicle_routes.items()):
#         color = colors[idx % len(colors)]
#         for (i, j) in route:
#             if i in nodes and j in nodes:
#                 net.add_edge(i, j, color=color, width=4, title=f"Vehicle {veh}")

#     # --- Draw transfer arcs ---
#     for (i, j) in transfer_paths:
#         if i in nodes and j in nodes:
#             net.add_edge(
#                 i, j,
#                 color="orange",
#                 width=3,
#                 dashes=True,
#                 title=f"Transfer {i}->{j}"
#             )

#     # --- Draw request arcs ---
#     for (i, j) in request_paths:
#         if i in nodes and j in nodes:
#             net.add_edge(
#                 i, j,
#                 color="rgba(255,255,255,0.6)",
#                 width=2,
#                 dashes=[5, 5],
#                 title=f"Request {i}->{j}"
#             )

#     # --- Save output ---
#     output_file = f"{model_name.replace(' ', '_').lower()}_routes.html"
#     net.write_html(output_file)
#     print(f"✅ Graph written to {output_file}")
