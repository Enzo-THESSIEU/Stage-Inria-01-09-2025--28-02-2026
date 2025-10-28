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
        tij = self.params['tij']
        C = self.sets['C']
        closest_node = None
        dist = float('inf')
        for transfer_node in C:
            if tij[self.base(node), self.base(transfer_node)] < dist:
                closest_node = transfer_node
                dist = tij[node, transfer_node]
        return closest_node
    
    def find_closest_transfer_dropoff(self, node):
        tij = self.params['tij']
        C = self.sets['C']
        closest_node = None
        dist = float('inf')
        for transfer_node in C:
            if tij[self.base(transfer_node), self.base(node)] < dist:
                closest_node = transfer_node
                dist = tij[transfer_node, node]
        return closest_node

    def define_PT_clusters(self):
        tij = self.params['tij']
        C = self.sets['C']
        P = self.sets['P']
        D = self.sets['D']
        preliminary_cluster = {}
        for p in P:
            closest_node = self.find_closest_transfer_pickup(p)
            preliminary_cluster[closest_node] += [p]
        for d in D:
            closest_node = self.find_closest_transfer_dropoff(d)
            preliminary_cluster[closest_node] += [d]
            


    ### === Change the routes that have been found to be understood by gurobi variables === ###
    
    def translate_to_gurobi_variables(self):
        x , v = self.vars_['x'], self.vars_['v']
        for vehicle in range(self.n_veh):
            for step in range(len(self.vehicle_route) - 1):
                x[vehicle, self.vehicle_route[step], self.vehicle_route[step + 1]] = 1


            
            
            








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


        
    
