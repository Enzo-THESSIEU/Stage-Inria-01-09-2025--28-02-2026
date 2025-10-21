import ast
import pandas as pd
import math
import plotly.graph_objects as go

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
                        if isinstance(elem[0], list) and len(elem[0]) == 2:
                            i, j = elem[0]
                        elif len(elem) >= 2 and isinstance(elem[0], int) and isinstance(elem[1], int):
                            i, j = elem[0], elem[1]
                        else:
                            continue
                        arcs.append((remap_node(i), remap_node(j)))
                veh_num = int(col.split()[1])
                vehicle_routes[veh_num] = arcs
            except Exception as e:
                print(f"⚠️ Error parsing {col}: {e}")
    return vehicle_routes


def display_vehicle_routes_from_idarp(csv_path = "C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes\\results_full_2h_1.csv"
                                      , model_name=None, title="IDARP Vehicle Routes (Plotly)"):
    """
    Visualize vehicle routes from an IDARP model CSV file using Plotly.
    """

    # --- Load CSV ---
    df = pd.read_csv(csv_path)

    # Select model
    if model_name:
        row = df[df["Model"] == model_name].iloc[0]
    else:
        row = df.iloc[0]
        model_name = row["Model"]

    print(f"🔍 Visualizing model: {model_name}")

    # --- Parse vehicle routes ---
    vehicle_routes = parse_vehicle_routes(row)

    # --- Create circular node layout ---
    nodes = {}
    for node_id in range(13):
        angle = 2 * math.pi * node_id / 13
        x = math.cos(angle)
        y = math.sin(angle)
        nodes[node_id] = (x, y)

    # --- Prepare Plotly traces ---
    fig = go.Figure()

    # Background faint arcs (optional)
    for i in nodes:
        for j in nodes:
            if i != j:
                x_coords = [nodes[i][0], nodes[j][0], None]
                y_coords = [nodes[i][1], nodes[j][1], None]
                fig.add_trace(go.Scatter(
                    x=x_coords, y=y_coords,
                    mode="lines",
                    line=dict(color="rgba(200,200,200,0.05)", width=1),
                    hoverinfo='none',
                    showlegend=False
                ))

    # --- Draw vehicle routes ---
    colors = ["#e74c3c", "#3498db", "#f1c40f", "#9b59b6", "#1abc9c", "#e67e22", "#2ecc71"]
    for idx, (veh, arcs) in enumerate(vehicle_routes.items()):
        color = colors[idx % len(colors)]
        for (i, j) in arcs:
            if i in nodes and j in nodes:
                x_coords = [nodes[i][0], nodes[j][0], None]
                y_coords = [nodes[i][1], nodes[j][1], None]
                fig.add_trace(go.Scatter(
                    x=x_coords, y=y_coords,
                    mode="lines+markers",
                    line=dict(color=color, width=3),
                    marker=dict(size=6, color=color),
                    name=f"Vehicle {veh}",
                    hovertext=[f"{i}→{j}"],
                    hoverinfo="text"
                ))

    # --- Add nodes as separate layer ---
    node_x = [coord[0] for coord in nodes.values()]
    node_y = [coord[1] for coord in nodes.values()]
    node_labels = list(nodes.keys())

    # Assign colors by group
    node_colors = []
    for n in node_labels:
        if n in [0, 9]:
            node_colors.append("#4CAF50")  # depots
        elif n in [1, 2, 3, 4]:
            node_colors.append("#e67e22")
        elif n in [5, 6, 7, 8]:
            node_colors.append("#9b59b6")
        else:
            node_colors.append("#3498db")

    fig.add_trace(go.Scatter(
        x=node_x, y=node_y,
        mode="markers+text",
        text=[f"{i}" for i in node_labels],
        textposition="top center",
        marker=dict(size=14, color=node_colors, line=dict(width=2, color='white')),
        name="Nodes"
    ))

    # --- Layout styling ---
    fig.update_layout(
        title=f"<b>{title}</b><br><sub>{model_name}</sub>",
        title_x=0.5,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(20,20,20,1)",
        showlegend=True,
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        height=750
    )

    # --- Save to HTML ---
    output_file = f"{model_name.replace(' ', '_').lower()}_routes_plotly.html"
    fig.write_html(output_file)
    print(f"✅ Plotly graph written to {output_file}")

    fig.show()


display_vehicle_routes_from_idarp()