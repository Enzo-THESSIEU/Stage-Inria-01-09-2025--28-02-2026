import numpy as np
import random as rd
import math as m
from collections import deque

class initial_solution:
    def __init__(self, m , vars_, sets, params, **IDARPoptions):
        self.m = m
        self. vars_ = vars_
        self.sets = sets
        self.params = params
        self.duplicate_transfers = IDARPoptions.get("duplicate_transfers", True)
        self.arc_elimination = IDARPoptions.get("arc_elimination", True)
        self.variable_substitution = IDARPoptions.get("variable_substitution", True)
        self.subtour_elimination = IDARPoptions.get("subtour_elimination", True)
        self.transfer_node_strengthening = IDARPoptions.get("transfer_node_strengthening", True)
        self.ev_constraints = IDARPoptions.get("ev_constraints", False)
        self.timetabled_departures = IDARPoptions.get("timetabled_departures", False)
        self.use_imjn = IDARPoptions.get("use_imjn", False)
        self.MoPS = IDARPoptions.get("MoPS", False)

        self.n_veh = len(self.sets['K'])

        # --- FIXED AND SAFE ---
        zero = self.sets['zeroDepot']

        self.vehicle_route = [[zero] for _ in range(self.n_veh)]
        self.vehicle_current_requests = [[] for _ in range(self.n_veh)]
        self.loc = [zero for _ in range(self.n_veh)]
        self.load = [0 for _ in range(self.n_veh)]

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
    
    def base_set(self, node_list):
        return set(self.base(i) for i in node_list)

    def order_by_earliest_time_window(self, unsorted_set):
        ei = self.params['ei']  # dict: {node: earliest_time}
        return sorted(unsorted_set, key=lambda i: ei[self.base(i)])

    def closest_dropoff(self, v):
        closest_drop = None
        min_time = float('inf')
        tij = self.params['tij']
        pair_pi_di = self.params['pair_pi_di']

        if len(self.vehicle_current_requests[v]) == 0:
            return None

        curr = self.base(self.loc[v])

        for req in self.vehicle_current_requests[v]:
            if self.use_imjn:
                drop_node = pair_pi_di[(req, 1)]
            else:
                drop_node = pair_pi_di[req]
            if drop_node is None:
                continue

            i = curr
            j = self.base(drop_node)

            # Key-safe lookup
            if (i, j) not in tij:
                continue

            if tij[i, j] < min_time:
                min_time = tij[i, j]
                closest_drop = drop_node

        return closest_drop

    def closest_vehicle(self, p):
        closest = None
        dist_closest = float('inf')
        tij = self.params['tij']

        for v in range(self.n_veh):
            i = self.base(self.loc[v])
            j = self.base(p)

            # Key-safe lookup
            if (i, j) not in tij:
                continue

            if tij[i, j] < dist_closest:
                closest = v
                dist_closest = tij[i, j]

        # Fallback if no arc found
        if closest is None:
            return 0

        return closest


    ### === Algorithm 2 in Overleaf: Closest Request Greedy Heuristic === ###
    def closest_request_greedy_heuristic(self):

        tij = self.params['tij']
        pair_pi_di = self.params['pair_pi_di']

        # ---- PHASE 1: Serve all pickups ----
        self.unserved_pickups = deque(self.order_by_earliest_time_window(self.sets['P']))
        print(f"Ordered pickups is {self.unserved_pickups}")

        while self.unserved_pickups:
            p = self.unserved_pickups.popleft()
            v = self.closest_vehicle(p)
            drop = self.closest_dropoff(v)
            # print(f"Pickup node is {p}, dropoff node is {drop}, closest vehicle to pickup is {v}")

            # Case 1: No dropoff available → always go pick up
            if drop is None:
                self.vehicle_route[v].append(p)
                self.vehicle_current_requests[v].append(self.base(p))
                self.loc[v] = p
                continue

            # Retrieve indices
            curr = self.base(self.loc[v])
            bp = self.base(p)
            bd = self.base(drop)

            # Safe lookup
            dist_p = tij[curr, bp] if (curr, bp) in tij else float('inf')
            dist_d = tij[curr, bd] if (curr, bd) in tij else float('inf')

            # Case 2: pickup is closer
            if dist_p < dist_d:
                # print(f"dist_p < dist_d: {dist_p} < {dist_d}")
                self.vehicle_route[v].append(p)
                self.vehicle_current_requests[v].append(self.base(p))
                self.loc[v] = p

            # Case 3: dropoff closer → serve dropoff and postpone pickup
            else:
                # print(f"dist_p > dist_d: {dist_p} > {dist_d}")
                self.vehicle_route[v].append(drop)
                self.loc[v] = drop
                self.unserved_pickups.appendleft(p)

                # remove its request id
                req = None
                for key, val in pair_pi_di.items():
                    if self.base(val) == bd:
                        req = self.base(key)
                        break
                if req in self.vehicle_current_requests[v]:
                    self.vehicle_current_requests[v].remove(req)

        # for v in range(self.n_veh):
            # print(f" For vehicle {v}, the current route is {self.vehicle_route[v]}, and the current requests onboard are {self.vehicle_current_requests[v]}")

        # ===================================================================
        # ---- PHASE 2: Serve all remaining dropoffs AFTER all pickups ------
        # ===================================================================

        for v in range(self.n_veh):

            while len(self.vehicle_current_requests[v]) > 0:

                drop = self.closest_dropoff(v)
                # print(f"The closest dropoff node for vehicle {v} is {drop}")
                if drop is None:
                    break  # shouldn't happen but safe guard

                bd = self.base(drop)
                self.vehicle_route[v].append(drop)
                self.loc[v] = drop

                # Remove dropped request
                req = None
                for key, val in pair_pi_di.items():
                    if self.base(val) == bd:
                        req = self.base(key)
                        break
                if req in self.vehicle_current_requests[v]:
                    self.vehicle_current_requests[v].remove(req)

        # ===================================================================
        # ---- PHASE 3: Send all vehicles to end depot -----------------------
        # ===================================================================

        end_depot = self.sets['endDepot']

        for v in range(self.n_veh):
            self.vehicle_route[v].append(end_depot)
            self.loc[v] = end_depot


    ### === Algorithm 3 in Overleaf: Earliest Time window Greedy heuristic for DARP === ###

    # def earliest_request_greedy_heuristic(self):
    #     tij = self.params['tij']
    #     self.unserved_P = self.order_by_earliest_time_window(self.sets['P'])
    #     self.unserved_D = self.order_by_earliest_time_window(self.sets['D'])
    #     self.unserved_P_and_D = self.unserved_P + self.unserved_D
    #     while self.unserved_P_and_D:
    #         next_node = self.unserved_P_and_D.pop()
    #         if next_node in self.unserved_D:
    #             pickup_node = next((k for k, v in self.params['pair_pi_di'].items() if v == next_node), None)
    #             vehicle_index = next((v for v, route in enumerate(self.vehicle_route) if pickup_node in route), None)
    #             self.vehicle_route[vehicle_index].append(next_node)
    #         else:
    #             next_vehicle = self.closest_vehicle(next_node)
    #             self.vehicle_route[next_vehicle].append(next_node)
    #             self.loc[next_vehicle] = next_node



    def earliest_request_greedy_heuristic(self):
        tij = self.params['tij']
        pair_pi_di = self.params['pair_pi_di']

        # Sort by earliest time window and convert to deque
        unserved_P = deque(self.order_by_earliest_time_window(self.sets['P']))
        unserved_D = deque(self.order_by_earliest_time_window(self.sets['D']))

        # Combine preserving order: P first, then D
        # (request dropoffs can be visited only after their pickups)
        self.unserved_P_and_D = deque(list(unserved_P) + list(unserved_D))

        while self.unserved_P_and_D:

            # Pop earliest node
            next_node = self.unserved_P_and_D.popleft()

            # CASE 1: next_node is a dropoff
            if next_node in unserved_D:

                # find pickup node for it
                pickup_node = next((k for k, v in pair_pi_di.items()
                                    if self.base(v) == self.base(next_node)), None)

                # find vehicle that served pickup
                vehicle_index = None
                for v_idx, route in enumerate(self.vehicle_route):
                    if pickup_node in route:
                        vehicle_index = v_idx
                        break

                if vehicle_index is None:
                    # no vehicle has served the pickup yet → cannot serve dropoff
                    # push it to the back of deque
                    self.unserved_P_and_D.append(next_node)
                    continue

                # append dropoff
                self.vehicle_route[vehicle_index].append(next_node)
                self.loc[vehicle_index] = next_node

                # remove from "carried" requests
                req = self.base(pickup_node)
                if req in self.vehicle_current_requests[vehicle_index]:
                    self.vehicle_current_requests[vehicle_index].remove(req)

            else:
            # CASE 2: next_node is a pickup

                # assign to closest vehicle
                next_vehicle = self.closest_vehicle(next_node)
                self.vehicle_route[next_vehicle].append(next_node)
                self.loc[next_vehicle] = next_node

                # mark request as being carried
                self.vehicle_current_requests[next_vehicle].append(self.base(next_node))


    ### === Algorithm 4 in Overleaf: Cluster-Based Intial Solution === ###

    def find_closest_transfer_pickup(self, node):
        C = self.base_set(self.sets['C'])
        tij = self.params['tij']
        closest_node = None
        dist = float('inf')
        for transfer_node in C:
            if tij[node, transfer_node] < dist:
                closest_node = transfer_node
                dist = tij[node, transfer_node]
        return closest_node
    
    def find_closest_transfer_dropoff(self, node):
        C = self.base_set(self.sets['C'])
        tij = self.params['tij']
        closest_node = None
        dist = float('inf')
        for transfer_node in C:
            if tij[transfer_node, node] < dist:
                closest_node = transfer_node
                dist = tij[transfer_node, node]
        return closest_node

    def define_PT_preliminary_clusters(self):
        P = self.base_set(self.sets['P'])
        D = self.base_set(self.sets['D'])
        if self.duplicate_transfers:
            C = sorted(self.params['Cr'][1])
        elif self.use_imjn:
            C = sorted(self.base_set(self.sets['C']))
        preliminary_clusters = {}
        for node in C:
            preliminary_clusters[node] = []
        for p in P:
            closest_node = self.find_closest_transfer_pickup(p)
            preliminary_clusters[closest_node] += [p]
        for d in D:
            closest_node = self.find_closest_transfer_dropoff(d)
            preliminary_clusters[closest_node] += [d]
        return preliminary_clusters

    def define_cluster_score(self, preliminary_clusters, C, tij):
        cluster_key, max_cluster = max(
                            preliminary_clusters.items(),
                            key=lambda item: len(item[1])
                        )
        
        max_tij = np.max(tij)

        C_score = {}
        if self.n_veh < len(C):
            for i_idx, i in enumerate(C):
                for j_idx, j in enumerate(C):
                    key = (i_idx, j_idx, i, j)
                    if key in C_score:
                        continue
                    if i_idx < j_idx:
                        C_score[key] = (len(preliminary_clusters[i]) + len(preliminary_clusters[j])) / len(max_cluster) + tij[i_idx][j_idx] / max_tij
                        
        return C_score
    
    def triangle_height(self, base, side_1, side_2):
        # a is the base
        s = (base + side_1 + side_2) / 2
        A = m.sqrt(s * (s - base) * (s - side_1) * (s - side_2))
        h = (2 * A) / base
        return h
    
    def calculate_distance_new_cluster_centre_to_node(self, tij, node_idx, cluster_i_idx, cluster_j_idx, w_i, w_j):
        dist_i_j = tij[cluster_i_idx, cluster_j_idx]
        dist_i_node = tij[cluster_i_idx, node_idx]
        dist_j_node = tij[cluster_j_idx, node_idx]
        dist_i_new_cluster = w_i * dist_i_j
        dist_j_new_cluster = w_j * dist_i_j
        triangle_height = self.triangle_height(dist_i_j, dist_i_node, dist_j_node)
        dist_cluster_middle_base = dist_i_j - dist_i_new_cluster
        dist_cluster_node = m.sqrt(dist_cluster_middle_base**2 + triangle_height**2)
        return dist_cluster_node

    def calculate_new_tij(self, link, tij, C, w_i, w_j):
        i_idx = link[0]
        j_idx = link[1]
        i = link[2]
        j = link[3]
        print("i_idx", i_idx)
        print("j_idx", j_idx)
        print("tij size", np.shape(tij))
        print("C is", C)
        new_dist = np.zeros(len(C))
        for k, node in enumerate(C):
            if k not in [i_idx, j_idx]:
                print("Index k", k)
                print("Node", node)
                new_dist[k] = self.calculate_distance_new_cluster_centre_to_node(tij, k, i_idx, j_idx, w_i, w_j)
                print(new_dist[k])
        new_dist = np.delete(new_dist, [i_idx, j_idx], axis=0)
        tij = np.delete(tij, [i_idx, j_idx], axis=0)
        tij = np.delete(tij, [i_idx, j_idx], axis=1)
        print("new dist column", new_dist)
        print("tij after deletion", tij)
        tij = np.hstack((tij, new_dist.reshape(-1, 1)))
        new_row = np.append(new_dist, 0)
        tij = np.vstack((tij, new_row))
        print("Final tij", tij)

        return tij

    def combine_clusters(self):
        preliminary_clusters = self.define_PT_preliminary_clusters()
        print(" ======================= \n")
        print("Preliminary Clsuters is ", preliminary_clusters)
        print(" ======================= \n")
        if self.duplicate_transfers:
            C = sorted(self.params['Cr'][1])
        elif self.use_imjn:
            C = sorted(self.base_set(self.sets['C']))
        tij = self.params['tij']

        # build dense tij submatrix for transfer nodes
        nC = len(C)
        tij_transfer_nodes = np.zeros((nC, nC))
        for a, i in enumerate(C):
            for b, j in enumerate(C):
                tij_transfer_nodes[a, b] = tij[self.base(i), self.base(j)]

        while len(C) > self.n_veh:
            C_score = self.define_cluster_score(preliminary_clusters, C, tij_transfer_nodes)
            link, link_val = min(C_score.items())
            w_i = len(preliminary_clusters[link[2]])
            w_j = len(preliminary_clusters[link[3]])
            preliminary_clusters[(link[2], link[3])] = preliminary_clusters[link[2]] + preliminary_clusters[link[3]]
            tij_transfer_nodes = self.calculate_new_tij(link, tij_transfer_nodes, C, w_i, w_j)
            C.remove(link[2])
            C.remove(link[3])
            C.append((link[2], link[3]))
            del C_score[link]
            del preliminary_clusters[link[2]]
            del preliminary_clusters[link[3]]
        return preliminary_clusters
    
    def find_nearest_transfer(self, node):
        C = self.sets["C"]
        # Assuming you have a cost/distance matrix `cij`
        return min(C, key=lambda t: self.params["tij"][node, t])
    
    def determine_opposite_node(self, node):
        pair_pi_di = self.params['pair_pi_di']
        P = self.sets['P']
        D = self.sets['D']
        if node in P:
            pickup_node = node
            dropoff_node = pair_pi_di[pickup_node]
        elif node in D:
            dropoff_node = node
            pickup_node = next(k for k, v in pair_pi_di.items() if v == dropoff_node)
        return pickup_node, dropoff_node
    
    def link_pickup_transfer_and_dropoff_transfer(self,node):
        preliminary_clusters = self.define_PT_preliminary_clusters()
        pickup_node, dropoff_node = self.determine_opposite_node(node)

        transfer_node_for_pickup_node = next(k for k, v in preliminary_clusters.items() if pickup_node in v)
        transfer_node_for_dropoff_node = next(k for k, v in preliminary_clusters.items() if dropoff_node in v)
        return transfer_node_for_pickup_node, transfer_node_for_dropoff_node
    
    def build_subproblem_for_each_cluster(self, cluster):
        print("Cluster is of format", cluster)
        pickup_nodes_for_cluster = []
        dropoff_nodes_for_cluster = []
        P = self.sets['P']
        D = self.sets['D']
        pair_pi_di = self.params['pair_pi_di']
        tij = self.params['tij']
        ei = self.params['ei']
        li = self.params['li']
        ei_cluster = []
        li_cluster = []
        len_cluster = len(cluster)
        for node in cluster:
            if node in P:
                pickup_node, dropoff_node = self.determine_opposite_node(node)
                transfer_node_for_pickup_node, transfer_node_for_dropoff_node = self.link_pickup_transfer_and_dropoff_transfer(node)
                pickup_nodes_for_cluster.append(node)
                dropoff_nodes_for_cluster.append(transfer_node_for_pickup_node)
                ei_cluster.append(ei[self.base(node)])
                li_cluster.append(li[self.base(dropoff_node)] - 
                                  tij[self.base(node), self.base(transfer_node_for_pickup_node)] - 
                                  tij[self.base(transfer_node_for_pickup_node), self.base(transfer_node_for_dropoff_node)] - 
                                  tij[self.base(transfer_node_for_dropoff_node), self.base(dropoff_node)])
            elif node in D:
                pickup_node, dropoff_node = self.determine_opposite_node(node)
                transfer_node_for_pickup_node, transfer_node_for_dropoff_node = self.link_pickup_transfer_and_dropoff_transfer(node)
                pickup_nodes_for_cluster.append(transfer_node_for_dropoff_node)
                dropoff_nodes_for_cluster.append(node)
                ei_cluster.append(ei[self.base(pickup_node)] + 
                                  tij[self.base(pickup_node), self.base(transfer_node_for_pickup_node)] + 
                                  tij[self.base(transfer_node_for_pickup_node), self.base(transfer_node_for_dropoff_node)])
                li_cluster.append(li[self.base(dropoff_node)] -  
                                  tij[self.base(transfer_node_for_dropoff_node), self.base(dropoff_node)])
        subproblem = [[] for _ in range(len_cluster)]   
        # Build array of arrays where the format is 
        # subproblem = [... [pickup node for subproblem, dropoff node for subproblem, ei time constraint for pickup node, id for subrequest] ...]
        for sub_prob_req_idx in range(len_cluster):
            subproblem[sub_prob_req_idx] = [pickup_nodes_for_cluster[sub_prob_req_idx], dropoff_nodes_for_cluster[sub_prob_req_idx], ei_cluster[sub_prob_req_idx], sub_prob_req_idx]
        subproblem = sorted(subproblem, key=lambda row: row[2])
        return deque(subproblem)
    
    def find_closest_dropoff_node_subproblem(self, pickups_in_vehicle, route, subproblem_copy):
        tij = self.params['tij']
        current_location = route[-1]
        low_dist = float('inf')
        closest_dropoff = pickups_in_vehicle[0]
        for pickup in pickups_in_vehicle:
            dropoff_node = next(
                                (item[1] for item in subproblem_copy if item[0] == pickup)
                                )
            if tij[self.base(current_location), self.base(dropoff_node)] < low_dist:
                low_dist = tij[self.base(current_location), self.base(dropoff_node)]
                closest_dropoff = dropoff_node
        return closest_dropoff
    
    def calculate_capacity_vehicle(self, pickups_in_vehicle):
        qr = self.params['qr']
        return sum(qr[pickup] for pickup in pickups_in_vehicle)
    
    def build_route_subproblem(self, subproblem):
        zeroDepot_node = self.sets['zeroDepot']
        Q = self.params['Q']
        qr = self.params['qr']
        tij = self.params['tij']
        route = [zeroDepot_node]
        pickups_in_vehicle = []
        subproblem_copy = subproblem.copy()
        while len(subproblem) > 0:
            next_node_list = subproblem.popleft()
            next_pickup = next_node_list[0]
            next_dropoff = None
            if len(pickups_in_vehicle) > 0:
                next_dropoff = self.find_closest_dropoff_node_subproblem(pickups_in_vehicle, route, subproblem_copy)
            if next_dropoff == None:   ### No request inside vehicle
                route.append(next_pickup)
                pickups_in_vehicle.append(next_pickup)
            elif self.calculate_capacity_vehicle(pickups_in_vehicle) + qr[next_pickup] > Q:  ### Will violate capacity constraint
                route.append(next_dropoff)
                subproblem.appendleft(next_node_list)
            elif tij[self.base(route[-1]), self.base(next_pickup)] < tij[self.base(route[-1]), self.base(next_dropoff)]:   ### pickup closer than dropoff
                route.append(next_pickup)
                pickups_in_vehicle.append(next_pickup)
            else:   ### Dropoff closer than pickup
                route.append(next_dropoff)
                subproblem.appendleft(next_node_list)
        return route
    
    def build_routes_for_clusters(self, clusters):
        routes = []
        for cluster_id, cluster in clusters.items():
            cluster_subproblem = self.build_subproblem_for_each_cluster(cluster)
            cluster_route = self.build_route_subproblem(cluster_subproblem)
            routes.append(cluster_route)
        return routes

    def Execute_Heuristic(self):
        clusters = self.combine_clusters()
        print(clusters)
        routes = self.build_routes_for_clusters(clusters)
        return routes

        
class translate_LNS_to_Gurobi:

    def __init__(self, vehicle_routes, m, vars_, sets, params, **options):
        self.vehicle_routes = vehicle_routes
        self.n_veh = len(self.vehicle_routes)
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params
        self.variable_substitution = options.get("variable_substitution")
        self.duplicate_transfers = options.get("duplicate_transfers")
        self.use_imjn = options.get("use_imjn")

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
  
    ### === Change the routes that have been found to be understood by gurobi variables === ###

    def intialise_all_gurobi_binary_variables_to_one(self):
        x, y, v, z, a = self.vars_['x'], self.vars_['y'], self.vars_['v'], self.vars_['z'], self.vars_['a']    
        for key in x:
            x[key].Start = 0

        for key in y:
            y[key].Start = 0

        for key in v:
            v[key].Start = 0

        for key in z:
            z[key].Start = 0

        for key in a:
            a[key].Start = 0

    def translate_to_gurobi_variables_x(self):
        """
        Takes self.vehicle_route[v] (a list of nodes for each vehicle)
        and assigns x or v variables = 1 accordingly.
        """

        print("Running translate_to_gurobi_variables_x")

        x = self.vars_.get("x", None)
        v = self.vars_.get("v", None)

        for k in range(self.n_veh):

            route = self.vehicle_routes[k]
            if len(route) <= 1:
                continue   # nothing to assign

            for idx in range(len(route) - 1):
                i = route[idx]
                j = route[idx + 1]

                if self.variable_substitution:
                    if (i, j) in v:
                        v[i, j].Start = 1
                        x[k, i, j].Start = 1
                        print(f"v[{i}, {j}] initialised to 1")
                        print(f"x[{k}, {i}, {j}] initialised to 1")
                else:
                    if (k, i, j) in x:
                        x[k, i, j].Start = 1
                        print(f"x[{k}, {i}, {j}] initialised to 1")

    # def translate_to_gurobi_variables_y_z_duplicate_transfers(self):
    #     """
    #     For each request r, reconstruct its full path including:
    #     - same-vehicle pickup → dropoff
    #     - or pickup → transfer → PT → transfer → dropoff

    #     Then set the corresponding Gurobi variables:
    #     - y[r,i,j] = 1 for vehicle arcs used by the request
    #     - z[r,t1,t2] = 1 for PT arcs between transfer nodes
    #     """

    #     y = self.vars_["y"]
    #     z = self.vars_["z"]

    #     pair_pi_di = self.params["pair_pi_di"]
    #     Cr = self.params["Cr"]  # Cr per request

    #     for r_pick, r_drop in pair_pi_di.items():

    #         r = r_pick  # request id = pickup node id

    #         # ---- STEP 1: find vehicle that visits the pickup ----
    #         v_pick = None
    #         idx_pick = None
    #         for v in range(self.n_veh):
    #             route = self.vehicle_routes[v]
    #             for idx, node in enumerate(route):
    #                 if self.base(node) == self.base(r_pick):
    #                     v_pick = v
    #                     idx_pick = idx
    #                     break
    #             if v_pick is not None:
    #                 break

    #         if v_pick is None:
    #             print(f"Request {r_pick} not served")
    #             continue  # request not served at all (should not happen)

    #         route1 = self.vehicle_routes[v_pick]

    #         # ---- STEP 2: search forward for the dropoff ----
    #         idx_drop = None
    #         for idx2 in range(idx_pick + 1, len(route1)):
    #             if self.base(route1[idx2]) == self.base(r_drop):
    #                 idx_drop = idx2
    #                 break

    #         # ---- CASE A: same vehicle serves pickup → dropoff ----
    #         if idx_drop is not None:
    #             for a in range(idx_pick, idx_drop):
    #                 i = self.base(route1[a])
    #                 j = self.base(route1[a+1])
    #                 if (r, i, j) in y:
    #                     y[r, i, j].Start = 1
    #             continue

    #         # ---- CASE B: request must transfer using PT ----
    #         # Find first transfer node in Cr on vehicle v_pick's route
    #         idx_t1 = None
    #         t1 = None
    #         for a in range(idx_pick + 1, len(route1)):
    #             node = self.base(route1[a])
    #             if node in Cr:
    #                 idx_t1 = a
    #                 t1 = node
    #                 break

    #         if t1 is None:
    #             print(f"Request {r_pick} has an incomplete route")
    #             print("The issue stems from pickup to transfer")
    #             continue  # should not happen if Cr is correct

    #         # Mark vehicle arcs from pickup → transfer t1
    #         for a in range(idx_pick, idx_t1):
    #             i = self.base(route1[a])
    #             j = self.base(route1[a+1])
    #             if (r, i, j) in y:
    #                 y[r, i, j].Start = 1

    #         # ---- STEP 3: find second vehicle that visits a Cr node ----
    #         v_drop = None
    #         idx_t2 = None
    #         t2 = None

    #         for v2 in range(self.n_veh):
    #             route2 = self.vehicle_route[v2]
    #             for a2, node in enumerate(route2):
    #                 node_b = self.base(node)
    #                 if node_b in Cr and node_b != t1:
    #                     v_drop = v2
    #                     idx_t2 = a2
    #                     t2 = node_b
    #                     break
    #             if v_drop is not None:
    #                 break

    #         if t2 is None:
    #             print(f"Request {r_pick} has an incomplete route")
    #             print("The issue stems from transfer to transfer")
    #             continue  # no second transfer → cannot serve

    #         # ---- STEP 4: PT arc between t1 → t2 ----
    #         if (t1, t2) in z:
    #             z[t1, t2].Start = 1

    #         # ---- STEP 5: continue vehicle v_drop from t2 → r_drop ----
    #         route2 = self.vehicle_route[v_drop]
    #         idx_final_drop = None

    #         for a2 in range(idx_t2 + 1, len(route2)):
    #             if self.base(route2[a2]) == self.base(r_drop):
    #                 idx_final_drop = a2
    #                 break

    #         if idx_final_drop is None:
    #             print(f"Request {r_pick} has an incomplete route")
    #             print("The issue stems from transfer to dropoff")
    #             continue  # incomplete second half

    #         for a2 in range(idx_t2, idx_final_drop):
    #             i = self.base(route2[a2])
    #             j = self.base(route2[a2+1])
    #             if (r, i, j) in y:
    #                 y[r, i, j].Start = 1

    # def translate_to_gurobi_variables_y_z_imjn(self):
    #     """
    #     Same as translate_to_gurobi_variables_y_z but adapted for IMJN node format,
    #     where transfer nodes are shared across requests.
        
    #     For a request that changes vehicles:
    #         - t1 = first transfer node in C visited AFTER pickup
    #         - t2 = last  transfer node in C visited BEFORE dropoff
    #     """

    #     y = self.vars_["y"]
    #     z = self.vars_["z"]
    #     pair_pi_di = self.params["pair_pi_di"]
    #     C = set(self.base(i) for i in self.sets["C"])

    #     for r_pick, r_drop in pair_pi_di.items():

    #         r = r_pick
    #         base_pick = self.base(r_pick)
    #         base_drop = self.base(r_drop)

    #         # === STEP 1 — find vehicle & position of pickup ===
    #         v_pick, idx_pick = None, None
    #         for v in range(self.n_veh):
    #             route = self.vehicle_routes[v]
    #             for idx, node in enumerate(route):
    #                 if self.base(node) == base_pick:
    #                     v_pick = v
    #                     idx_pick = idx
    #                     break
    #             if v_pick is not None:
    #                 break
    #         if v_pick is None:
    #             print(f"Request {self.base(r_pick)} not served")
    #             continue  # request not served

    #         route1 = self.vehicle_routes[v_pick]

    #         # === STEP 2 — does same vehicle visit dropoff? ===
    #         idx_drop_same = None
    #         for a in range(idx_pick + 1, len(route1)):
    #             if self.base(route1[a]) == base_drop:
    #                 idx_drop_same = a
    #                 break

    #         if idx_drop_same is not None:
    #             # SAME VEHICLE: mark arcs pickup → dropoff
    #             for a in range(idx_pick, idx_drop_same):
    #                 i = self.base(route1[a])
    #                 j = self.base(route1[a+1])
    #                 if (r, i, j) in y:
    #                     y[r, i, j].Start = 1
    #             continue

    #         # === STEP 3 — Request must transfer ===
    #         # t1 = first transfer node after pickup
    #         t1, idx_t1 = None, None
    #         for a in range(idx_pick + 1, len(route1)):
    #             node_b = self.base(route1[a])
    #             if node_b in C:
    #                 t1 = node_b
    #                 idx_t1 = a
    #                 break
    #         if t1 is None:
    #             print(f"Request {self.base(r_pick)} has an incomplete route")
    #             print("The issue stems from pickup to transfer")
    #             continue  # no transfer found → invalid route

    #         # mark pickup → t1
    #         for a in range(idx_pick, idx_t1):
    #             i = self.base(route1[a])
    #             j = self.base(route1[a+1])
    #             if (r, i, j) in y:
    #                 y[r, i, j].Start = 1

    #         # === STEP 4 — find t2: last transfer before dropoff ===
    #         v_drop, idx_drop, t2, idx_t2 = None, None, None, None

    #         for v2 in range(self.n_veh):
    #             route2 = self.vehicle_route[v2]
    #             # first find dropoff on this vehicle
    #             for a2, node in enumerate(route2):
    #                 if self.base(node) == base_drop:
    #                     v_drop = v2
    #                     idx_drop = a2
    #                     break
    #             if v_drop is not None:
    #                 # now search BACKWARDS from idx_drop for last transfer
    #                 for a2 in range(idx_drop - 1, -1, -1):
    #                     node_b = self.base(route2[a2])
    #                     if node_b in C and node_b != t1:
    #                         t2 = node_b
    #                         idx_t2 = a2
    #                         break
    #                 break

    #         if t2 is None:
    #             print(f"Request {self.base(r_pick)} has an incomplete route")
    #             print("The issue stems from transfer to transfer")
    #             continue  # no second transfer → cannot serve request

    #         # === STEP 5 — assign z[r,t1,t2] ===
    #         # if (r, t1, t2) in z:
    #         #     z[r, t1, t2].Start = 1

    #         # === STEP 6 — from t2 → dropoff on v_drop ===
    #         route2 = self.vehicle_route[v_drop]
    #         for a2 in range(idx_t2, idx_drop):
    #             i = self.base(route2[a2])
    #             j = self.base(route2[a2+1])
    #             if (r, i, j) in y:
    #                 y[r, i, j].Start = 1

    def translate_to_gurobi_variables_y_z_duplicate_transfers(self):
        """
        For each request r, reconstruct its full path including:
        - direct:  pickup -> dropoff on the same vehicle
        - via PT: pickup -> t1 (on vehicle v_pick),
                z[t1,t2] (PT arc),
                t2 -> dropoff (on vehicle v_drop)

        This is fully *route-based*, so it works for:
        - heuristics that only generate direct routes
        - heuristics that only generate PT routes
        - heuristics that mix both, per request.
        """

        y = self.vars_["y"]
        z = self.vars_.get("z", {})

        pair_pi_di = self.params["pair_pi_di"]      # {pickup: dropoff}
        Cr = self.params["Cr"]                  # either dict r->iterable(C_r) or global iterable

        for r_pick, r_drop in pair_pi_di.items():
            r = r_pick
            base_pick = self.base(r_pick)
            base_drop = self.base(r_drop)

            # Per-request transfer set Cr(r)
            if isinstance(Cr, dict):
                Cr_r = {self.base(c) for c in Cr[r]}
            else:
                Cr_r = {self.base(c) for c in Cr}

            # === STEP 1 — find vehicle & position of pickup ===
            v_pick, idx_pick = None, None
            for v in range(self.n_veh):
                route = self.vehicle_routes[v]
                for idx, node in enumerate(route):
                    if self.base(node) == base_pick:
                        v_pick = v
                        idx_pick = idx
                        break
                if v_pick is not None:
                    break

            if v_pick is None:
                print(f"[y/z dup] Request {base_pick} not served (no pickup)")
                continue

            route1 = self.vehicle_routes[v_pick]

            # === STEP 2 — does the same vehicle also visit the dropoff? ===
            idx_drop_same = None
            for a in range(idx_pick + 1, len(route1)):
                if self.base(route1[a]) == base_drop:
                    idx_drop_same = a
                    break

            if idx_drop_same is not None:
                # DIRECT CASE: same vehicle pickup -> dropoff
                for a in range(idx_pick, idx_drop_same):
                    i = self.base(route1[a])
                    j = self.base(route1[a + 1])
                    if (r, i, j) in y:
                        y[r, i, j].Start = 1
                        print(f"y[{r}, {i}, {j}] initialised to 1")
                continue  # go to next request

            # === STEP 3 — via transfer: find t1 (first Cr node after pickup on v_pick) ===
            t1, idx_t1 = None, None
            for a in range(idx_pick + 1, len(route1)):
                node_b = self.base(route1[a])
                if node_b in Cr_r:
                    t1 = node_b
                    idx_t1 = a
                    break

            if t1 is None:
                print(f"[y/z dup] Request {base_pick} incomplete: pickup->transfer (no t1 found)")
                continue

            # Mark pickup -> t1 arcs on v_pick
            for a in range(idx_pick, idx_t1):
                i = self.base(route1[a])
                j = self.base(route1[a + 1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1
                    print(f"y[{r}, {i}, {j}] initialised to 1")

            # === STEP 4 — find t2 and v_drop: last Cr node before dropoff, on some vehicle ===
            v_drop, idx_t2, t2, idx_drop = None, None, None, None

            for v2 in range(self.n_veh):
                route2 = self.vehicle_routes[v2]

                # positions of dropoff on this route
                drop_positions = [idx for idx, node in enumerate(route2)
                                if self.base(node) == base_drop]
                if not drop_positions:
                    continue
                # take first dropoff visit on this vehicle
                local_idx_drop = drop_positions[0]

                # candidate transfer nodes before that dropoff, excluding t1
                transfer_positions = [idx for idx, node in enumerate(route2)
                                    if idx < local_idx_drop
                                    and self.base(node) in Cr_r
                                    and self.base(node) != t1]
                if not transfer_positions:
                    continue

                # t2 = last transfer before dropoff
                local_idx_t2 = max(transfer_positions)
                candidate_t2 = self.base(route2[local_idx_t2])

                v_drop = v2
                idx_t2 = local_idx_t2
                t2 = candidate_t2
                idx_drop = local_idx_drop
                break  # assume unique good match

            if t2 is None or v_drop is None:
                print(f"[y/z dup] Request {base_pick} incomplete: transfer->transfer (no t2 found)")
                continue

            # === STEP 5 — PT arc between t1 -> t2 ===
            # z can be indexed as (t1,t2) or (r,t1,t2) depending on your model.
            # Here we first try (r,t1,t2), then (t1,t2).
            if (r, t1, t2) in z:
                z[r, t1, t2].Start = 1
            elif (t1, t2) in z:
                z[t1, t2].Start = 1

            # === STEP 6 — mark arcs t2 -> dropoff on v_drop ===
            route2 = self.vehicle_routes[v_drop]
            for a2 in range(idx_t2, idx_drop):
                i = self.base(route2[a2])
                j = self.base(route2[a2 + 1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1
                    print(f"y[{r}, {i}, {j}] initialised to 1")

    def translate_to_gurobi_variables_y_z_imjn(self):
        """
        IMJN version (shared transfer nodes):
        - auto-detect direct vs via-transfer per request
        - if via-transfer:
                t1 = first transfer in C after pickup on v_pick
                t2 = last  transfer in C before dropoff on v_drop
        - only sets y; z is intentionally ignored for IMJN.
        """

        print("Running correctly translate_to_gurobi_variables_y_z_imjn")

        y = self.vars_["y"]
        pair_pi_di = self.params["pair_pi_di"]
        C = self.sets["C"]   # global set of transfer nodes

        for r_pick, r_drop in pair_pi_di.items():
            request = self.base(r_pick)
            base_pick = self.base(r_pick)
            base_drop = self.base(r_drop)

            # === STEP 1 — find vehicle & position of pickup ===
            v_pick, idx_pick = None, None
            for v in range(self.n_veh):
                route = self.vehicle_routes[v]
                for idx, node in enumerate(route):
                    if self.base(node) == base_pick:
                        v_pick = v
                        idx_pick = idx
                        break
                if v_pick is not None:
                    break

            if v_pick is None:
                print(f"[y/z imjn] Request {base_pick} not served (no pickup)")
                continue

            route1 = self.vehicle_routes[v_pick]

            # === STEP 2 — direct case: same vehicle visits dropoff ===
            idx_drop_same = None
            for a in range(idx_pick + 1, len(route1)):
                if self.base(route1[a]) == base_drop:
                    idx_drop_same = a
                    break

            if idx_drop_same is not None:
                # print(f"Direct route found for request {request}")
                # DIRECT: pickup -> dropoff on v_pick
                for a in range(idx_pick, idx_drop_same):
                    i = route1[a]
                    j = route1[a + 1]
                    # print((i,j))
                    if (request, i, j) in y:
                        y[request, i, j].Start = 1
                        print(f"y[{request}, {i}, {j}] initialised to 1")
                    else:
                        print(f"{(request, i, j)} not in y")
                continue  # next request

            # === STEP 3 — via-transfer: t1 = first C node after pickup on v_pick ===
            t1, idx_t1 = None, None
            for a in range(idx_pick + 1, len(route1)):
                node_b = route1[a]
                if node_b in C:
                    t1 = node_b
                    idx_t1 = a
                    break

            if t1 is None:
                print(f"[y/z imjn] Request {base_pick} incomplete: pickup->transfer (no t1 found)")
                continue

            # mark arcs pickup -> t1 on v_pick
            for a in range(idx_pick, idx_t1):
                i = route1[a]
                j = route1[a + 1]
                if (request, i, j) in y:
                    y[request, i, j].Start = 1
                    print(f"y[{request}, {i}, {j}] initialised to 1")
                else:
                    print(f"{(request, i, j)} not in y")

            # === STEP 4 — find v_drop, t2 and idx_drop ===
            v_drop, idx_drop, t2, idx_t2 = None, None, None, None

            for v2 in range(self.n_veh):
                route2 = self.vehicle_routes[v2]

                # find first dropoff position on this vehicle
                idx_drop_candidate = None
                for a2, node in enumerate(route2):
                    if self.base(node) == base_drop:
                        idx_drop_candidate = a2
                        break
                if idx_drop_candidate is None:
                    continue

                # scan backwards to get last transfer before that dropoff, different from t1
                for a2 in range(idx_drop_candidate - 1, -1, -1):
                    node_b = route2[a2]
                    if node_b in C and node_b != t1:
                        v_drop = v2
                        idx_drop = idx_drop_candidate
                        t2 = node_b
                        idx_t2 = a2
                        break

                if v_drop is not None:
                    break

            if t2 is None or v_drop is None:
                print(f"[y/z imjn] Request {base_pick} incomplete: transfer->dropoff (no t2 found)")
                continue

            # === STEP 5 — mark arcs t2 -> dropoff on v_drop ===
            route2 = self.vehicle_routes[v_drop]
            for a2 in range(idx_t2, idx_drop):
                i = route2[a2]
                j = route2[a2 + 1]
                if (request, i, j) in y:
                    y[request, i, j].Start = 1
                    print(f"y[{request}, {i}, {j}] initialised to 1")
                else:
                    print(f"{(request, i, j)} not in y")



class destroy_operators:
    def __init__(self, m, vars_, sets, params, **options):
        self.m = m
        self. vars_ = vars_
        self.sets = sets
        self.params = params

        # Extract nested option groups
        IDARPoptions = options.get("IDARPoptions", {})
        Heuristic_options = options.get("Heuristic_options", {})

        # Model options
        self.duplicate_transfers = IDARPoptions.get("duplicate_transfers", True)
        self.arc_elimination = IDARPoptions.get("arc_elimination", True)
        self.variable_substitution = IDARPoptions.get("variable_substitution", True)
        self.subtour_elimination = IDARPoptions.get("subtour_elimination", True)
        self.transfer_node_strengthening = IDARPoptions.get("transfer_node_strengthening", True)
        self.ev_constraints = IDARPoptions.get("ev_constraints", False)
        self.timetabled_departures = IDARPoptions.get("timetabled_departures", False)
        self.use_imjn = IDARPoptions.get("use_imjn", False)
        self.MoPS = IDARPoptions.get("MoPS", False)

        # LNS options
        initial_solution_options = Heuristic_options.get("initial_solution_options", {})
        removal_heuristic_options = Heuristic_options.get("removal_heuristic_options", {})

        # Removal Heuistic Choice Options
        random_removal_options = removal_heuristic_options.get("random_removal_heuristic_options", {})
        worst_removal_options = removal_heuristic_options.get("worst_removal_heuristic_options", {})
        related_removal_options = removal_heuristic_options.get("related_removal_heuristic_options", {})

        # Random Removal options
        self.random_removal_rate = random_removal_options.get("random_removal_rate", 0.2)

        # Worst Removal Options
        self.worst_removal_rate = worst_removal_options.get("worst_removal_options", 0.2)

        # Related Removal Options
        self.remove_n_closest = related_removal_options.get("remove_n_closest", False)
        self.n_closest_number = related_removal_options.get("n_closest_number", 1)
        self.remove_closest_by_r = related_removal_options.get("remove_closest_by_r", False)
        self.radius_r = related_removal_options.get("radius_r", 1)
        
        # Link Model -> Heuristics
        self.routes = self.build_routes_from_solution()
    
    def build_routes_from_solution(self):
        """
        Builds vehicle routes directly from solved Gurobi variables x[k,i,j].
        Returns:
            routes: dict {vehicle k: [route sequence including depot(s)]}
        """
        x = self.vars_['x']             # Gurobi var dict: x[k,i,j]
        K = self.sets['K']              # set/list of vehicles
        N = self.sets['N']              # all nodes
        A = self.sets['A']              # all arcs (i,j)
        zeroDepot = self.sets['zeroDepot']
        endDepot  = self.sets['endDepot']

        routes = [[] for k in range(len(K))]

        for k in K:
            # find arcs used by vehicle k
            arcs_k = [(i, j) for (kk, i, j) in x.keys() if kk == k and x[kk, i, j].X > 0.5]

            if not arcs_k:
                continue  # vehicle unused

            # build route by chaining arcs
            route = [zeroDepot]
            current_node = route[-1]

            while True:
                next_nodes = [j for (i, j) in arcs_k if i == current_node]
                if not next_nodes:
                    break
                next_node = next_nodes[0]
                route.append(next_node)
                current_node = next_node
                if current_node == endDepot:
                    break

            routes[k] = route

        return routes
    
    ### === Base Functions === ###

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
    
    def base_set(self, node_list):
        return set(self.base(i) for i in node_list)
    
    def aggregate_pickups_dropoffs(self, removed_pickups):
        removed_dropoffs = [self.params['pair_pi_di'][p] for p in removed_pickups]
        return removed_pickups + removed_dropoffs
    
    def remove_list_routes(self, removal_list):
        # for each vehicle route
        for r_idx, route in enumerate(self.routes):
            # build a new route excluding removed nodes
            new_route = [node for node in route if node not in removal_list]
            self.routes[r_idx] = new_route
    
    ### === Random Removal Heuristic === ###

    def random_removal_choice(self):
        P = list(self.sets['P'])
        n_remove = max(1, int(len(P) * self.random_removal_rate))  # remove at least 1
        
        removed_pickups = rd.sample(P, n_remove)
        remaining_pickups = [p for p in P if p not in removed_pickups]
        
        return removed_pickups, remaining_pickups
    
    def random_removal_route(self):
        # randomly choose which pickups to remove
        removed_pickups = self.random_removal_choice()
        removal_list = self.aggregate_pickups_dropoffs(removed_pickups)
        
        self.remove_list_routes(removal_list)

        return removal_list
 
    ### === Worst Removal Heuristic === ###

    def compute_objective_value(self, routes, w=[1, 0, 0, 0]):
        """
        Computes the objective value numerically from route lists (no Gurobi vars).
        Used during LNS to evaluate solution quality.

        w = [w1, w2, w3, w4]  (same weighting convention as in the model)
        """
        cij = self.params["cij"]
        tij = self.params["tij"]
        di  = self.params["di"]
        pair_pi_di = self.params["pair_pi_di"]

        # === f1: total travel cost (sum of distances along all arcs) ===
        f1 = 0
        for route in routes:
            for i, j in zip(route[:-1], route[1:]):
                f1 += cij[self.base(i), self.base(j)]

        # === f2 and f3: time window and service penalties (simplified) ===
        # For LNS evaluation, we can approximate using total delay between pickups & dropoffs
        f2 = 0
        f3 = 0
        for p, d in pair_pi_di.items():
            # find which route contains this pair
            for route in routes:
                if p in route and d in route:
                    idx_p = route.index(p)
                    idx_d = route.index(d)
                    # approximate time between pickup and dropoff as sum of tij along that segment
                    travel_time = sum(tij[self.base(route[r]), self.base(route[r + 1])] for r in range(idx_p, idx_d))
                    service_time = di[self.base(p)]
                    # direct shortest path time for normalization
                    direct_time = tij[self.base(p), self.base(d)]
                    f2 += (travel_time - service_time - direct_time)
                    f3 += (travel_time - service_time) / (direct_time + 1e-6)  # avoid div by zero
                    break  # skip to next pair

        # === f4: optional MoPS term (if applicable) ===
        if self.MoPS:
            D_M = self.sets['D_M']
            f4 = sum(
                1 for route in routes for i, j in zip(route[:-1], route[1:])
                if i in D_M
            )
        else:
            f4 = 0

        # === Final weighted sum ===
        total_obj = w[0] * f1 + w[1] * f2 + w[2] * f3 - w[3] * f4
        return total_obj

    def remove_request(self, request_node):
        removal_nodes = [request_node, self.params['pair_pi_di'][request_node]]
        routes = self.routes
        for r_idx, route in enumerate(self.routes):
            # build a new route excluding removed nodes
            new_route = [node for node in route if node not in removal_nodes]
            routes[r_idx] = new_route
        return routes
    
    def worst_removal(self):
        P = self.sets['P']
        objective_values= {}
        for p in P:
            new_routes = self.remove_request(p)
            obj_val = self.compute_objective_value(new_routes)
            objective_values[p] = obj_val
        total_requests = len(P)
        requests_removed = 0
        requests_to_remove = []
        while requests_removed < max(1, total_requests * self.random_removal_rate):
            key, val = max(objective_values.items())
            del objective_values[key]
            requests_removed += 1
            requests_to_remove.append(key)
        nodes_to_remove = self.aggregate_pickups_dropoffs(requests_to_remove)
        self.remove_list_routes(nodes_to_remove)

    ### === Related Removal Heuristic === ###

    def find_n_closest_neighbours(self, n, node):
        zeroDepot, endDepot = self.base(self.sets['zeroDepot']), self.base(self.sets['endDepot'])
        tij_node = self.params['tij'][self.base(node)]

        # Make sure that depot nodes aren't removed from the routes
        tij_node[zeroDepot] = float('inf')
        tij_node[endDepot] = float('inf')

        n_closest = sorted(enumerate(tij_node), key=lambda x: x[1])[:n]
        return n_closest

    def find_closest_neighbours_by_r(self, r, node):
        zeroDepot, endDepot = self.base(self.sets['zeroDepot']), self.base(self.sets['endDepot'])
        tij_node = self.params['tij'][self.base(node)]
        r_dist = max(tij_node) * r

        # Make sure that depot nodes aren't removed from the routes
        tij_node[zeroDepot] = float('inf')
        tij_node[endDepot] = float('inf')

        closest_by_r = []
        for i in range(tij_node):
            if tij_node[i] < r_dist:
                closest_by_r.append(i)

        return closest_by_r
    
    def remove_closest_neighbours(self):
        zeroDepot, endDepot = self.base(self.sets['zeroDepot']), self.base(self.sets['endDepot'])
        potential_nodes = [n for n in self.sets['N'] if n not in [zeroDepot, endDepot]]
        node = rd.random(potential_nodes)
        if self.remove_n_closest:
            removal_list = self.find_n_closest_neighbours(self.n_closest_number, node)
        elif self.remove_closest_by_r:
            removal_list = self.find_closest_neighbours_by_r(self.radius_r, node)
        self.remove_list_routes(removal_list)










        
            








# class LNS_IDARP:
#     def __init__(self, m , vars_, sets, params, **IDARPoptions, **LNSoptions):
#         self.m = m
#         self. vars_ = vars_
#         self.sets = sets
#         self.params = params
#         self.duplicate_transfers = IDARPoptions.get("duplicate_transfers", True)
#         self.arc_elimination = IDARPoptions.get("arc_elimination", True)
#         self.variable_substitution = IDARPoptions.get("variable_substitution", True)
#         self.subtour_elimination = IDARPoptions.get("subtour_elimination", True)
#         self.transfer_node_strengthening = IDARPoptions.get("transfer_node_strengthening", True)
#         self.ev_constraints = IDARPoptions.get("ev_constraints", False)
#         self.timetabled_departures = IDARPoptions.get("timetabled_departures", False)
#         self.use_imjn = IDARPoptions.get("use_imjn", False)
#         self.MoPS = IDARPoptions.get("MoPS", False)
#         self.iter = LNSoptions.get("iterations", 1)
    
#     def do_neighbourhood_search(self, s):
#         total_iter = 0
#         s1 = s
#         while total_iter <= self.iter:
#             s1 = self.destroy_heuristic(self, s1)


        
    
