import jsonpickle
jsonpickle.backend.remove_backend('simplejson')
jsonpickle.backend.remove_backend('rapidjson')
jsonpickle.backend.remove_backend('ujson')

import json
import ast
import re
from pathlib import Path
from pyvis.network import Network


class HallPlotter:

    # ================================================================
    # INIT
    # ================================================================
    def __init__(self, logs_dir="logs"):
        self.logs_dir = Path(logs_dir)
        if not self.logs_dir.exists():
            raise FileNotFoundError(f"Logs directory not found: {logs_dir}")

        self.log_path = self._get_latest_log_file()
        self.data = self._extract_last_json_object()


    # ================================================================
    # GET NEWEST LOG
    # ================================================================
    def _get_latest_log_file(self):
        files = sorted(self.logs_dir.glob("*.log"), reverse=True)
        if not files:
            raise FileNotFoundError("No .log files found in logs/")
        return files[0]


    # ================================================================
    # EXTRACT LAST JSON BLOCK
    # ================================================================
    def _extract_last_json_object(self):
        """
        Extracts the final JSON object from inside the .log file.
        """
        txt = self.log_path.read_text()

        # regex capturing nested { ... } JSON blocks
        json_blocks = re.findall(r"\{(?:[^{}]|(?:\{[^{}]*\}))*\}", txt, flags=re.DOTALL)

        if not json_blocks:
            raise ValueError("No JSON objects found in log.")

        last_json = json_blocks[-1]

        try:
            parsed = json.loads(last_json)
        except Exception as e:
            raise ValueError(f"Failed parsing JSON:\n{e}\nJSON was:\n{last_json}")

        return parsed


    # ================================================================
    # COLLAPSE IMJN NODE → BASE NODE
    # ================================================================
    def _to_base_node(self, imjn):
        """
        Convert (i,m) node to its base i.
        Used to collapse transfer/charging layers.
        """
        if not isinstance(imjn, tuple) or len(imjn) != 2:
            return None
        return imjn[0]


    # ================================================================
    # PARSE VEHICLE ROUTES
    # ================================================================
    def _parse_vehicle_routes(self):
        """
        Convert "Veh1": "[[( (i,m),(j,n) ), ...], ...]" into
        a list of base arcs (i -> j), colouring transfers/chargers separately.
        """
        routes = {}

        for key, raw in self.data.items():
            if not key.startswith("Veh"):
                continue
            if not isinstance(raw, str):
                continue

            # Safely parse Python list from string
            try:
                steps = ast.literal_eval(raw)
            except:
                continue

            arcs = []
            types = []   # keep node types to identify charging nodes

            for elem in steps:
                # elem must be: [[(i,m),(j,n)], node, type, req, ...]
                if not isinstance(elem, list) or len(elem) < 3:
                    continue

                arc = elem[0]
                node_type = elem[2]

                if not (isinstance(arc, list) and len(arc) == 2):
                    continue

                prev_node = arc[0]   # (i,m)
                next_node = arc[1]   # (j,n)

                i = self._to_base_node(prev_node)
                j = self._to_base_node(next_node)

                if i is None or j is None:
                    continue

                arcs.append((i, j))
                types.append((i, node_type))  # store node type for special colouring

            routes[key] = {
                "arcs": arcs,
                "types": types
            }

        return routes


    # ================================================================
    # FIXED HALL LAYOUT (TRANSFER + CHARGING CO-LOCATED)
    # ================================================================
    def _hall_layout(self):
        """
        Charge stations (13,14,15) are plotted right next to transfers (10,11,12).
        """

        layout = {}

        # Depots
        layout[0] = (0, +10)
        layout[9] = (0, -10)

        # Pickups
        layout[1] = (+10, +20)
        layout[2] = (-40, +250)
        layout[3] = (-40, +200)
        layout[4] = (-200, 0)

        # Dropoffs
        layout[5] = (-80, +250)
        layout[6] = (-150, -50)
        layout[7] = (-150, 0)
        layout[8] = (-80, +200)

        # Transfer nodes
        layout[10] = (-60, +200)
        layout[11] = (-10, -20)
        layout[12] = (-200, -50)

        # Charging nodes: shifted slightly so they appear visually adjacent
        layout[13] = (-60 + 10, +200)   # near transfer 10
        layout[14] = (-10 + 10, -20)    # near transfer 11
        layout[15] = (-200 + 10, -50)   # near transfer 12

        return layout


    # ================================================================
    # COLOR SCHEME FOR NODES
    # ================================================================
    def _color_node(self, n):
        if n in [0, 9]:
            return "#4CAF50"  # depots
        if 1 <= n <= 4:
            return "#E67E22"  # pickups
        if 5 <= n <= 8:
            return "#9B59B6"  # dropoffs
        if n in [10, 11, 12]:
            return "#3498DB"  # transfers
        if n in [13, 14, 15]:
            return "#1ABC9C"  # charging stations
        return "#7F8C8D"


    # ================================================================
    # MAIN PLOTTING FUNCTION
    # ================================================================
    def plot(self):
        model_name = self.data.get("Model", "HallPlot")
        routes = self._parse_vehicle_routes()
        layout = self._hall_layout()

        # PyVis network
        net = Network(
            height="900px",
            width="100%",
            directed=True,
            bgcolor="white",
            font_color="black"
        )
        net.toggle_physics(False)

        # ---- Add nodes ----
        for n, (x, y) in layout.items():
            net.add_node(
                n,
                label=str(n),
                x=x,
                y=y,
                shape="dot",
                size=12,
                physics=False,
                color=self._color_node(n)
            )

        # ---- Vehicle paths ----
        colors = ["#E74C3C", "#2980B9"]  # Veh1 red, Veh2 blue

        for idx, (veh, content) in enumerate(routes.items()):
            c = colors[idx % len(colors)]
            arcs = content["arcs"]

            for (i, j) in arcs:
                if i in layout and j in layout:
                    net.add_edge(
                        i, j,
                        color=c,
                        width=4,
                        arrows="to",
                        title=f"{veh}: {i}→{j}",
                    )

        # Save output
        out = f"hall_plot_{model_name}.html"
        net.write_html(out)

        print(f"\n✅ Plot saved as {out}")
        print(f"📄 Using log: {self.log_path}")


# ================================================================
# RUN WHEN DIRECTLY CALLED
# ================================================================
if __name__ == "__main__":
    HallPlotter().plot()


