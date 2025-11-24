import json
import math
from pathlib import Path
from pyvis.network import Network
import ast


class HallPlotter:

    def __init__(self, logs_dir="logs"):
        self.logs_dir = Path(logs_dir)
        if not self.logs_dir.exists():
            raise FileNotFoundError(f"Logs directory not found: {logs_dir}")

        self.log_path = self._get_latest_log_file()
        self.data = self._extract_last_json_object()

    # ----------------------------------------------------------
    # 1. FIND NEWEST .log FILE
    # ----------------------------------------------------------
    def _get_latest_log_file(self):
        log_files = sorted(self.logs_dir.glob("*.log"), reverse=True)
        if not log_files:
            raise FileNotFoundError("No .log files found in logs/")
        return log_files[0]

    # ----------------------------------------------------------
    # 2. EXTRACT LAST JSON OBJECT INSIDE .log
    # ----------------------------------------------------------
    def _extract_last_json_object(self):
        """
        The .log contains lines like:
        2025-11-20 11:36:13 - INFO - { ...json... }

        This function extracts ALL JSON blocks { ... } and returns the last.
        """

        text = self.log_path.read_text()

        # Extract all { ... } blocks using a REGEX that captures balanced JSON objects
        json_blocks = re.findall(r"\{(?:[^{}]|(?:\{[^{}]*\}))*\}", text, flags=re.DOTALL)

        if not json_blocks:
            raise ValueError("No JSON objects found inside log file.")

        last_json_str = json_blocks[-1]

        try:
            parsed = json.loads(last_json_str)
        except json.JSONDecodeError as e:
            raise ValueError(f"Failed to parse final JSON object.\nError: {e}\nJSON was:\n{last_json_str}")

        return parsed

    # def _load_last_entry(self):
    #     """Load the LAST logged experiment entry."""
    #     with open(self.json_path, "r") as f:
    #         payload = json.load(f)

    #     # Each entry is {"title": "...", "data": {...}}
    #     last_entry = payload[-1]["data"]
    #     return last_entry

    # --------------------------------------------------------
    # DATA PARSING
    # --------------------------------------------------------
    def _normalize_node(self, n):
        DUP_MAP = {
            13: 10, 16: 10, 19: 10,
            14: 11, 17: 11, 20: 11,
            15: 12, 18: 12, 21: 12
        }
        return DUP_MAP.get(n, n)

    def _parse_vehicle_routes(self):
        routes = {}

        for key, raw in self.data.items():
            print("key", key)
            print("raw", raw)

            if not ("Veh" in key or "Vehicle" in key):
                continue

            if not isinstance(raw, str):
                continue

            try:
                parsed = ast.literal_eval(raw)
                print(parsed)
            except Exception as e:
                print(f"Failed to parse vehicle string for {key}: {e}")
                continue

            arcs = []

            # Each element in vehicle list = 1 route step
            for elem in parsed:
                # Must be a list of length >=1
                if not isinstance(elem, list) or len(elem) == 0:
                    continue

                # elem[0] MUST be: [(prev_node,m), (next_node,n)]
                if not (isinstance(elem[0], list) and len(elem[0]) == 2):
                    continue

                prev_node = elem[0][0]  # (i,k)
                next_node = elem[0][1]  # (j,k)

                # sanity check
                if not (isinstance(prev_node, tuple) and isinstance(next_node, tuple)):
                    continue

                i = self._normalize_node(prev_node[0])
                j = self._normalize_node(next_node[0])

                arcs.append((i, j))

            routes[key] = arcs

        return routes



    # --------------------------------------------------------
    # FIXED HALL LAYOUT (MATCHING FIGURE 2)
    # --------------------------------------------------------
    def _hall_layout(self):
        nodes = {}
        # --- DEPOTS (top & bottom center) ---
        nodes[0] = (0, +10)   # Origin depot
        nodes[9] = (0, -10)   # Destination depot

        # --- PICK-UPS (left side cluster) ---
        nodes[1] = (+10, +20)
        nodes[2] = (-40, +250)
        nodes[3] = (-40, +200)
        nodes[4] = (-200, 0)

        # --- DROPOFFS (right side cluster) ---
        nodes[5] = (-80, +250)
        nodes[6] = (-150, -50)
        nodes[7] =  (-150, 0)
        nodes[8] = (-80, +200)

        # --- TRANSFER NODES (center-right triangle) ---
        nodes[10] = (-60, +200)
        nodes[11] = (-10, -20)
        nodes[12] = (-200, -50)
        return nodes

    # --------------------------------------------------------
    # COLOR SCHEME
    # --------------------------------------------------------
    def _color_node(self, n):
        if n in [0, 9]:
            return "#4CAF50"  # Green depots
        if 1 <= n <= 4:
            return "#E67E22"  # Orange pickups
        if 5 <= n <= 8:
            return "#9B59B6"  # Purple dropoffs
        if n in [10, 11, 12]:
            return "#3498DB"  # Blue transfer nodes
        return "#95A5A6"      # Fallback grey

    # --------------------------------------------------------
    # MAIN PLOT FUNCTION
    # --------------------------------------------------------
    def plot(self):
        model_name = self.data.get("Model", "HallPlot")
        vehicle_routes = self._parse_vehicle_routes()

        print(vehicle_routes)

        layout = self._hall_layout()

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
                x=x, y=y,
                shape="dot",
                size=10,
                physics=False,
                color=self._color_node(n)
            )

        # ---- Vehicle Routes ----
        colors = ["#E74C3C", "#2980B9", "#F1C40F", "#2ECC71"]

        for idx, (veh, arcs) in enumerate(vehicle_routes.items()):
            c = colors[idx % len(colors)]
            for (i, j) in arcs:
                if i in layout and j in layout:
                    net.add_edge(
                        i, j,
                        color=c,
                        width=4,
                        title=f"{veh}: {i}->{j}"
                    )

        out = f"hall_plot_{model_name}.html"
        net.write_html(out)

        print(f"✅ Plot saved as: {out}")
        print(f"📄 Using JSON log: {self.json_path}")


# --------------------------------------------------------
# Optional: Direct execution
# --------------------------------------------------------
if __name__ == "__main__":
    HallPlotter().plot()

