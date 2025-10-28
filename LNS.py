import numpy as np

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
                    if i_idx < j_idx:
                        C_score[i, j] = (len(preliminary_clusters[i]) + len(preliminary_clusters[j])) / max_cluster + tij[i, j] / max_tij
        return C_score
    
    def calculate_new_tij(link):
        i = link[0]
        j = link[1]
        tij_1 = tij[:]


    
    def combine_clusters(self):
        preliminary_clusters = self.define_PT_preliminary_clusters()
        C = self.base_set(self.sets['C'])
        tij_transfer_nodes = self.params['tij'][min(C):max(C), min(C):max(C)]
        while len(C) > self.n_veh:
            C_score = self.define_cluster_score(preliminary_clusters, C, tij_transfer_nodes)
            link, link_val = min(C_score.items())
            del C_score[link]
            del preliminary_clusters[link[0]]
            del preliminary_clusters[link[1]]
            C.remove(link[0], link[1])
            C.add(link)





        


    ### === Change the routes that have been found to be understood by gurobi variables === ###
    
    def translate_to_gurobi_variables(self):
        x , v = self.vars_['x'], self.vars_['v']
        for vehicle in range(self.n_veh):
            for step in range(len(self.vehicle_route) - 1):
                x[vehicle, self.vehicle_route[step], self.vehicle_route[step + 1]] = 1


            
class destroy_operators:
    def __init__(self, m, vars_, sets, params, **IDARPoptions):
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

    def base(self, i):
        return i[0] if isinstance(i, tuple) else i
    
    def base_set(self, node_list):
        return set(self.base(i) for i in node_list)

    def random_removal(self, lambda_val):
        num_removals = int(lambda_val * len(self.params['P']))
        for a in range(num_removals):

        
            








class LNS_IDARP:
    def __init__(self, m , vars_, sets, params, **IDARPoptions, **LNSoptions):
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
        self.iter = LNSoptions.get("iterations", 1)
    
    def do_neighbourhood_search(self, s):
        total_iter = 0
        s1 = s
        while total_iter <= self.iter:
            s1 = self.destroy_heuristic(self, s1)


        
    
