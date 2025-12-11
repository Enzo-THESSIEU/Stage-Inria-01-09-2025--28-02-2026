import numpy as np


class data_extractor:
    def __init__(self, file_path, mrt_factor):
        self.file = file_path
        self.mrt_factor = mrt_factor

    def extract_node_info(self):
        nodes = []

        with open(self.file, "r") as f:
            # --- Read first line: n_vehicles n_requests n_terminals ---
            first = f.readline().strip().split()
            self.n_vehicles = int(first[0])
            self.n_requests = int(first[1])
            self.n_transfers = int(first[2])

            # --- Compute index ranges like in the C++ code ---
            pickup_max = self.n_requests
            dropoff_max = 2 * self.n_requests
            depot_max = 2 * self.n_requests + self.n_vehicles
            transfer_max = 2 * self.n_requests + self.n_vehicles + self.n_transfers

            # --- Read all remaining lines ---
            for line in f:
                parts = line.strip().split()
                if len(parts) < 8:
                    continue  # skip empty / invalid lines

                node_id = int(parts[0])
                x = float(parts[1])       # ignored but readable
                y = float(parts[2])       # ignored but readable
                lw = float(parts[3])
                uw = float(parts[4])
                mobile = int(parts[5])    # ignored here
                wheelchair = int(parts[6])# ignored here
                service_time = float(parts[7])

                # --- Detect node type ---
                if node_id < pickup_max:
                    node_type = "pickup"
                    request = node_id

                elif node_id < dropoff_max:
                    node_type = "dropoff"
                    request = node_id - self.n_requests

                elif node_id < depot_max:
                    node_type = "depot"
                    request = None

                elif node_id < transfer_max:
                    node_type = "transfer"
                    request = None

                else:
                    # Should not happen
                    node_type = "unknown"
                    request = None

                capacity = 1

                nodes.append([
                    node_id, node_type, request,
                    lw, uw, service_time, capacity
                ])

        return nodes

    def extract_travel_matrix(self, pt_speed):

        with open(self.file, "r") as f:
            self.num_nodes = 2 * self.n_requests + self.n_vehicles + self.n_transfers

            tij = np.zeros((self.num_nodes, self.num_nodes))
            x = np.zeros(self.num_nodes)
            y = np.zeros(self.num_nodes)

            for line in f:
                parts = line.strip().split()
                if len(parts) < 8:
                    continue  # skip empty / invalid lines

                idx = int(parts[0])
                x[idx] = float(parts[1])
                y[idx] = float(parts[2])

            for i_idx in range(self.num_nodes):
                for j_idx in range(self.num_nodes):

                    # ---------------------------------------------------------
                    # CASE 1 — Both nodes are PT terminals
                    # ---------------------------------------------------------
                    if i_idx >= 2 * self.n_requests + self.n_vehicles and j_idx >= 2 * self.n_requests + self.n_vehicles:

                        # Coordinates
                        xi, yi = x[i_idx], y[i_idx]
                        xj, yj = x[j_idx], y[j_idx]

                        # -----------------------------------------------------
                        # PT EXCLUSION RULES (illegal shortcuts)
                        # Mirrors of the C++ code:
                        # -----------------------------------------------------
                        illegal_x_sym = (xi == -xj) and (abs(yi) != 10) and (abs(yj) != 10)
                        illegal_y_sym = (yi == -yj) and (abs(xi) != 10) and (abs(xj) != 10)

                        if i_idx == j_idx:
                            tij[i_idx, j_idx] = 0.0

                        elif illegal_x_sym or illegal_y_sym:
                            # Forbidden PT arc
                            tij[i_idx, j_idx] = 10000.0

                        else:
                            # Manhattan PT travel time × speed factor
                            tij[i_idx, j_idx] = (abs(xi - xj) + abs(yi - yj)) * pt_speed

                    # ---------------------------------------------------------
                    # CASE 2 — At least one node is NOT a PT terminal
                    # Standard Euclidean distance
                    # ---------------------------------------------------------
                    else:
                        tij[i_idx, j_idx] = np.sqrt((x[j_idx] - x[i_idx])**2 +
                                                    (y[j_idx] - y[i_idx])**2)

            
        return tij
    
    def compute_maximum_ride_time(self, tij):
        mrt = np.zeros(self.n_requests)
        for i in range(self.n_requests):
            mrt[i] = tij[i, self.n_requests+i] * self.mrt_factor
        return mrt

    def tighten_time_windows(self, nodes, tij):
        mrt = self.compute_maximum_ride_time(tij)
        for node in nodes:
            if node[1] == "dropoff":
                node_id = node[0]
                request = node[2]
                lw = node[3]
                uw = node[4]
                lw_pickup = lw - mrt[request]
                uw_pickup = uw - tij[request, node_id]
                nodes[node_id - self.n_requests][3] = lw_pickup
                nodes[node_id - self.n_requests][4] = uw_pickup
        return nodes




                








