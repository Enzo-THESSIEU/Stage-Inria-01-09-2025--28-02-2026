import pandas as pd
import ast
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import Graph_Builder_3 as g

# --- Load results from CSV ---
results = {}
df = pd.read_csv("results.csv")
for _, row in df.iterrows():
    name = row["Model"]
    objval = float(row["Objective"])
    bound = float(row["Lower Bound"])
    v1_route = ast.literal_eval(row["Vehicle 1 Route"])
    v2_route = ast.literal_eval(row["Vehicle 2 Route"])
    cpu_time = float(row["CPU Time (s)"])
    elapsed = float(row["Elapsed Time (s)"])
    results[name] = (objval, bound, v1_route, v2_route, cpu_time, elapsed)

tij = g.data['tij_original']

# Build master graph layout
G_master = nx.DiGraph()
for _, (_, _, v1_route, v2_route, _, _) in results.items():
    for route in [v1_route, v2_route]:
        for step in route:
            u, v = step[0]
            G_master.add_edge(u, v)

pos = nx.spring_layout(G_master, seed=42)

# Helper for axis suffix
def axis_suffix(idx):
    return "" if idx == 1 else str(idx)

# Grid: Veh1 | Veh2 | Obj/LB | CPU/Elapsed | Name
fig = make_subplots(
    rows=len(results), cols=5,
    column_widths=[0.25, 0.25, 0.15, 0.2, 0.2],
    horizontal_spacing=0.05,
    subplot_titles=[]
)

plot_idx = 0
row_idx = 0
for name, (objval, bound, v1_route, v2_route, cpu_time, elapsed) in results.items():
    row_idx += 1

    for vidx, route in enumerate([v1_route, v2_route], start=1):
        plot_idx += 1
        col = vidx
        color = "red" if vidx == 1 else "blue"

        G = nx.DiGraph()
        for step in route:
            u, v = step[0]
            w = tij[(u, v)]
            G.add_edge(u, v, weight=w)

        sub_pos = {n: pos[n] for step in route for n in step[0]}

        # Edge lines + arrows
        edge_x, edge_y = [], []
        for u, v in G.edges():
            x0, y0 = sub_pos[u]
            x1, y1 = sub_pos[v]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]

            subplot_idx = (row_idx - 1) * 5 + col   # 5 columns in your grid
            xref = "x" if subplot_idx == 1 else f"x{subplot_idx}"
            yref = "y" if subplot_idx == 1 else f"y{subplot_idx}"

            fig.add_annotation(
                x=x1, y=y1, ax=x0, ay=y0,
                xref=xref, yref=yref,
                axref=xref, ayref=yref,
                arrowhead=3, arrowsize=1.2, arrowwidth=1.5,
                arrowcolor=color,
                row=row_idx, col=col
            )

        fig.add_trace(
            go.Scatter(x=edge_x, y=edge_y,
                       line=dict(width=2, color=color),
                       mode='lines', hoverinfo='none'),
            row=row_idx, col=col
        )

        # Nodes
        node_x, node_y, node_text, node_hover = [], [], [], []
        for step in route:
            node = step[0][0]
            service_time = step[4]
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(str(node))
            node_hover.append(f"Node {node}<br>Service time: {service_time}")

        fig.add_trace(
            go.Scatter(x=node_x, y=node_y,
                       mode='markers+text',
                       text=node_text,
                       textposition="bottom center",
                       marker=dict(size=18, color='white',
                                   line=dict(width=2, color='black')),
                       hovertext=node_hover, hoverinfo='text'),
            row=row_idx, col=col
        )
        i = 0
        for node in route:
            fig.add_trace(
                go.Scatter(
                    x=[0], y=[-i], mode="text",
                    text=[f"Node: {node[1]}<\r>"
                          f"Time: {node[4]}"],
                    hoverinfo="skip"
                ),
                row=row_idx, col=3+vidx
            )
            i+=1

    # ===== Metrics as 3 columns =====
    # # Column 4: Vehicle 1
    # for node in v1_route:
    #     fig.add_trace(
    #         go.Scatter(
    #             x=[0], y=[0], mode="text",
    #             text=[f"Node number: {node[1]}"],
    #             text=[f"Service Time: {node[4]}"]
    #             hoverinfo="skip"
    #         ),
    #         row=row_idx, col=3
    #     )

    # # Column 5: Vehicle 2
    # fig.add_trace(
    #     go.Scatter(
    #         x=[0], y=[0], mode="text",
    #         text=[f"CPU: {cpu_time:.2f}s<br>Elapsed: {elapsed:.2f}s"],
    #         hoverinfo="skip"
    #     ),
    #     row=row_idx, col=4
    # )

    # Column 3: Model Name and Metrics
    fig.add_trace(
        go.Scatter(
            x=[0], y=[0], mode="text",
            text=[f"<b>{name}</b><br>"
                f"Obj: {objval:.1f}<br>"
                f"LB: {bound:.1f}<br>"
                f"CPU: {cpu_time:.2f}s<br>"
                f"Elapsed: {elapsed:.2f}s"],
            hoverinfo="skip"
        ),
        row=row_idx, col=3
    )

# Clean up layout
fig.update_layout(
    showlegend=False,
    height=300*len(results), width=1800,
    margin=dict(t=40, b=40),
    plot_bgcolor="light gray",
    # paper_bgcolor="blue"
)

# Make all axes invisible (so only arrows/nodes/text remain)
# for axis in fig.layout:
#     if axis.startswith("xaxis") or axis.startswith("yaxis"):
#         fig.layout[axis].visible = False
fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False)
fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False)

# fig.write_image("plot.png") 
fig.show()
