import numpy as np
import random as rd
import math as m

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
        self.vehicle_route = np.array([ self.params['zeroDepot_node'] ], self.n_veh)
        self.vehicle_current_requests = np.array(self.n_veh, dtype=np.array)
        self.loc = np.array([ self.params['zeroDepot_node'] ], self.n_veh)
        self.load = np.zeros(self.n_veh)

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
    
    def base_set(self, node_list):
        return set(self.base(i) for i in node_list)

    def order_by_earliest_time_window(self, unsorted_set):
        ei = self.params['ei']  # dict: {node: earliest_time}
        return sorted(unsorted_set, key=lambda i: ei[self.base(i)])

    def closest_dropoff(self, v):
        closest_drop = None
        min_time = float('inf')  # Initialise with infinity to ensure any valid time will be smaller
        tij = self.params['tij']

        for req in self.vehicle_current_requests[v]:
            drop_node = self.params['pair_pi_di'][req]
            current_time = tij[self.base(self.loc[v]), self.base(drop_node)]

            if current_time < min_time:
                min_time = current_time
                closest_drop = drop_node

        return closest_drop
    
    def closest_vehicle(self, p):
        closest = None
        dist_closest = float('inf')
        tij = self.params['tij']
        for v in range(self.n_veh):
            if tij[self.base(self.loc[v]) , self.base(p)] < dist_closest:
                closest = v
        return closest
    
    ### === Algorithm 2 in Overleaf: Closest Request Greedy Heuristic === ###

    def closest_request_greedy_heuristic(self):
        tij = self.params['tij']
        self.unserved_pickups = self.order_by_earliest_time_window(self.sets['P'])
        while self.unserved_pickups:
            p = self.unserved_pickups.pop()
            closest_vehicle = self.closest_vehicle(p)
            closest_dropoff = self.closest_dropoff(closest_vehicle)
            if tij[self.base(self.loc[closest_vehicle]), self.base(p)] < tij[self.base(self.loc[closest_vehicle]), self.base(closest_dropoff)]:
                self.loc[closest_vehicle] = p
                self.vehicle_route[closest_vehicle].append(p)
            else:
                self.loc[closest_vehicle] = closest_dropoff
                self.vehicle_route[closest_vehicle].append(closest_dropoff)
                self.unserved_pickups.add(p)

    ### === Algorithm 3 in Overleaf: Earliest Time window Greedy heuristic for DARP === ###

    def earliest_request_greedy_heuristic(self):
        tij = self.params['tij']
        self.unserved_P = self.order_by_earliest_time_window(self.sets['P'])
        self.unserved_D = self.order_by_earliest_time_window(self.sets['D'])
        self.unserved_P_and_D = self.unserved_P + self.unserved_D
        while self.unserved_P_and_D:
            next_node = self.unserved_P_and_D.pop()
            if next_node in self.unserved_D:
                pickup_node = next((k for k, v in self.params['pair_pi_di'].items() if v == next_node), None)
                vehicle_index = next((v for v, route in enumerate(self.vehicle_route) if pickup_node in route), None)
                self.vehicle_route[vehicle_index].append(next_node)
            else:
                next_vehicle = self.closest_vehicle(next_node)
                self.vehicle_route[next_vehicle].append(next_node)
                self.loc[next_vehicle] = next_node

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
        preliminary_clusters = {}
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
        max_tij = max(tij[min(C):max(C), min(C):max(C)])
        C_score = {}
        if self.n_veh < len(C):
            for i, i_idx in enumerate(C):
                for j, j_idx in enumerate(C):
                    if C_score[i, j, i_idx, j_idx]: continue
                    else:
                        if i_idx < j_idx:
                            C_score[i, j, i_idx, j_idx] = (len(preliminary_clusters[i]) + len(preliminary_clusters[j])) / max_cluster + tij[i, j] / max_tij
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
        i = link[0]
        j = link[1]
        i_idx = link[2]
        j_idx = link[3]
        tij = np.delete(tij, [i_idx, j_idx], axis=0)
        tij = np.delete(tij, [i_idx, j_idx], axis=1)
        new_dist = np.zeros(len(C))
        for k in len(C):
            new_dist[k] = self.calculate_distance_new_cluster_centre_to_node(tij, k, i_idx, j_idx, w_i, w_j)
        tij = np.concatenate((tij, new_dist), axis = 0)
        tij = np.concatenate((tij, new_dist + [0]), axis = 1)
        return tij

    def combine_clusters(self):
        preliminary_clusters = self.define_PT_preliminary_clusters()
        C = self.base_set(self.sets['C'])
        tij_transfer_nodes = self.params['tij'][min(C):max(C), min(C):max(C)]
        while len(C) > self.n_veh:
            C_score = self.define_cluster_score(preliminary_clusters, C, tij_transfer_nodes)
            link, link_val = min(C_score.items())
            w_i = len(preliminary_clusters[link[0]])
            w_j = len(preliminary_clusters[link[1]])
            C.remove(link[0], link[1])
            C.add(link)
            preliminary_clusters[link] = preliminary_clusters[link[0]] + preliminary_clusters[link[1]]
            tij_transfer_nodes = self.calculate_new_tij(link, tij_transfer_nodes, C, w_i, w_j)
            del C_score[link]
            del preliminary_clusters[link[0]]
            del preliminary_clusters[link[1]]
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
        subproblem = np.array(len_cluster)   
        # Build array of arrays where the format is 
        # subproblem = [... [pickup node for subproblem, dropoff node for subproblem, ei time constraint for pickup node, id for subrequest] ...]
        for sub_prob_req_idx in range(len_cluster):
            subproblem[sub_prob_req_idx] = np.array(pickup_nodes_for_cluster[sub_prob_req_idx], dropoff_nodes_for_cluster[sub_prob_req_idx], ei_cluster[sub_prob_req_idx], sub_prob_req_idx)
        return subproblem
    
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
            min_idx = np.argmin([row[2] for row in subproblem])
            next_node_list = subproblem.pop(min_idx)
            next_pickup = next_node_list[0]
            next_dropoff = None
            if len(pickups_in_vehicle) > 0:
                next_dropoff = self.find_closest_dropoff_node_subproblem(pickups_in_vehicle, route, subproblem_copy)
            if next_dropoff == None:   ### No request inside vehicle
                route.append(next_pickup)
                pickups_in_vehicle.append(next_pickup)
            elif self.calculate_capacity_vehicle + qr[next_pickup] > Q:  ### Will violate capacity constraint
                route.append(next_dropoff)
                subproblem.add(next_node_list)
            elif tij[self.base(route[-1]), self.base(next_pickup)] < tij[self.base(route[-1]), self.base(next_dropoff)]:   ### pickup closer than dropoff
                route.append(next_pickup)
                pickups_in_vehicle.append(next_pickup)
            else:   ### Dropoff closer than pickup
                route.append(next_dropoff)
                subproblem.add(next_node_list)
        return route
    
    def build_routes_for_clusters(self, clusters):
        routes = []
        for cluster in clusters:
            cluster_subproblem = self.build_subproblem_for_each_cluster(cluster)
            cluster_route = self.build_route_subproblem(cluster_subproblem)
            routes.append(cluster_route)
        return routes


        
class translate_LNS_to_Gurobi:

    def __init__(self, vehicle_routes):
        self.vehicle_routes = vehicle_routes
        self.n_veh = len(self.vehicle_routes)

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
  
    ### === Change the routes that have been found to be understood by gurobi variables === ###
    
    def translate_to_gurobi_variables_x(self):
        """
        Takes self.vehicle_route[v] (a list of nodes for each vehicle)
        and assigns x or v variables = 1 accordingly.
        """

        x = self.vars_.get("x", None)
        v = self.vars_.get("v", None)

        for k in range(self.n_veh):

            route = self.vehicle_route[k]
            if len(route) <= 1:
                continue   # nothing to assign

            for idx in range(len(route) - 1):
                i = self.base(route[idx])
                j = self.base(route[idx + 1])

                # Standard case: x[k,i,j] are used
                if not self.variable_substitution:
                    try:
                        self.m.addConstr(x[k, i, j] == 1)
                    except KeyError:
                        pass  # arc not in model → ignore silently

                # Variable substitution: v[i,j] used instead of x
                else:
                    try:
                        self.m.addConstr(v[i, j] == 1)
                    except KeyError:
                        pass

    def translate_to_gurobi_variables_y_z_duplicate_transfers(self):
        """
        For each request r, reconstruct its full path including:
        - same-vehicle pickup → dropoff
        - or pickup → transfer → PT → transfer → dropoff

        Then set the corresponding Gurobi variables:
        - y[r,i,j] = 1 for vehicle arcs used by the request
        - z[r,t1,t2] = 1 for PT arcs between transfer nodes
        """

        y = self.vars_["y"]
        z = self.vars_["z"]

        pair_pi_di = self.params["pair_pi_di"]
        Cr = self.params["Cr"]  # Cr per request

        for r_pick, r_drop in pair_pi_di.items():

            r = r_pick  # request id = pickup node id

            # ---- STEP 1: find vehicle that visits the pickup ----
            v_pick = None
            idx_pick = None
            for v in range(self.n_veh):
                route = self.vehicle_routes[v]
                for idx, node in enumerate(route):
                    if self.base(node) == self.base(r_pick):
                        v_pick = v
                        idx_pick = idx
                        break
                if v_pick is not None:
                    break

            if v_pick is None:
                print(f"Request {r_pick} not served")
                continue  # request not served at all (should not happen)

            route1 = self.vehicle_routes[v_pick]

            # ---- STEP 2: search forward for the dropoff ----
            idx_drop = None
            for idx2 in range(idx_pick + 1, len(route1)):
                if self.base(route1[idx2]) == self.base(r_drop):
                    idx_drop = idx2
                    break

            # ---- CASE A: same vehicle serves pickup → dropoff ----
            if idx_drop is not None:
                for a in range(idx_pick, idx_drop):
                    i = self.base(route1[a])
                    j = self.base(route1[a+1])
                    if (r, i, j) in y:
                        y[r, i, j].Start = 1
                continue

            # ---- CASE B: request must transfer using PT ----
            # Find first transfer node in Cr on vehicle v_pick's route
            idx_t1 = None
            t1 = None
            for a in range(idx_pick + 1, len(route1)):
                node = self.base(route1[a])
                if node in Cr:
                    idx_t1 = a
                    t1 = node
                    break

            if t1 is None:
                print(f"Request {r_pick} has an incomplete route")
                print("The issue stems from pickup to transfer")
                continue  # should not happen if Cr is correct

            # Mark vehicle arcs from pickup → transfer t1
            for a in range(idx_pick, idx_t1):
                i = self.base(route1[a])
                j = self.base(route1[a+1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1

            # ---- STEP 3: find second vehicle that visits a Cr node ----
            v_drop = None
            idx_t2 = None
            t2 = None

            for v2 in range(self.n_veh):
                route2 = self.vehicle_route[v2]
                for a2, node in enumerate(route2):
                    node_b = self.base(node)
                    if node_b in Cr and node_b != t1:
                        v_drop = v2
                        idx_t2 = a2
                        t2 = node_b
                        break
                if v_drop is not None:
                    break

            if t2 is None:
                print(f"Request {r_pick} has an incomplete route")
                print("The issue stems from transfer to transfer")
                continue  # no second transfer → cannot serve

            # ---- STEP 4: PT arc between t1 → t2 ----
            if (t1, t2) in z:
                z[t1, t2].Start = 1

            # ---- STEP 5: continue vehicle v_drop from t2 → r_drop ----
            route2 = self.vehicle_route[v_drop]
            idx_final_drop = None

            for a2 in range(idx_t2 + 1, len(route2)):
                if self.base(route2[a2]) == self.base(r_drop):
                    idx_final_drop = a2
                    break

            if idx_final_drop is None:
                print(f"Request {r_pick} has an incomplete route")
                print("The issue stems from transfer to dropoff")
                continue  # incomplete second half

            for a2 in range(idx_t2, idx_final_drop):
                i = self.base(route2[a2])
                j = self.base(route2[a2+1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1

    def translate_to_gurobi_variables_y_z_imjn(self):
        """
        Same as translate_to_gurobi_variables_y_z but adapted for IMJN node format,
        where transfer nodes are shared across requests.
        
        For a request that changes vehicles:
            - t1 = first transfer node in C visited AFTER pickup
            - t2 = last  transfer node in C visited BEFORE dropoff
        """

        y = self.vars_["y"]
        z = self.vars_["z"]
        pair_pi_di = self.params["pair_pi_di"]
        C = set(self.base(i) for i in self.sets["C"])

        for r_pick, r_drop in pair_pi_di.items():

            r = r_pick
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
                print(f"Request {self.base(r_pick)} not served")
                continue  # request not served

            route1 = self.vehicle_routes[v_pick]

            # === STEP 2 — does same vehicle visit dropoff? ===
            idx_drop_same = None
            for a in range(idx_pick + 1, len(route1)):
                if self.base(route1[a]) == base_drop:
                    idx_drop_same = a
                    break

            if idx_drop_same is not None:
                # SAME VEHICLE: mark arcs pickup → dropoff
                for a in range(idx_pick, idx_drop_same):
                    i = self.base(route1[a])
                    j = self.base(route1[a+1])
                    if (r, i, j) in y:
                        y[r, i, j].Start = 1
                continue

            # === STEP 3 — Request must transfer ===
            # t1 = first transfer node after pickup
            t1, idx_t1 = None, None
            for a in range(idx_pick + 1, len(route1)):
                node_b = self.base(route1[a])
                if node_b in C:
                    t1 = node_b
                    idx_t1 = a
                    break
            if t1 is None:
                print(f"Request {self.base(r_pick)} has an incomplete route")
                print("The issue stems from pickup to transfer")
                continue  # no transfer found → invalid route

            # mark pickup → t1
            for a in range(idx_pick, idx_t1):
                i = self.base(route1[a])
                j = self.base(route1[a+1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1

            # === STEP 4 — find t2: last transfer before dropoff ===
            v_drop, idx_drop, t2, idx_t2 = None, None, None, None

            for v2 in range(self.n_veh):
                route2 = self.vehicle_route[v2]
                # first find dropoff on this vehicle
                for a2, node in enumerate(route2):
                    if self.base(node) == base_drop:
                        v_drop = v2
                        idx_drop = a2
                        break
                if v_drop is not None:
                    # now search BACKWARDS from idx_drop for last transfer
                    for a2 in range(idx_drop - 1, -1, -1):
                        node_b = self.base(route2[a2])
                        if node_b in C and node_b != t1:
                            t2 = node_b
                            idx_t2 = a2
                            break
                    break

            if t2 is None:
                print(f"Request {self.base(r_pick)} has an incomplete route")
                print("The issue stems from transfer to transfer")
                continue  # no second transfer → cannot serve request

            # === STEP 5 — assign z[r,t1,t2] ===
            # if (r, t1, t2) in z:
            #     z[r, t1, t2].Start = 1

            # === STEP 6 — from t2 → dropoff on v_drop ===
            route2 = self.vehicle_route[v_drop]
            for a2 in range(idx_t2, idx_drop):
                i = self.base(route2[a2])
                j = self.base(route2[a2+1])
                if (r, i, j) in y:
                    y[r, i, j].Start = 1


            
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


        
    
