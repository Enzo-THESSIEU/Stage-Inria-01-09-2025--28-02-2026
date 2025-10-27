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
        self.vehicle_route = np.array(self.n_veh, dtype=np.array)
        self.vehicle_current_requests = np.array(self.n_veh, dtype=np.array)
        self.vehicle_closest_dropoff = np.array(self.n_veh, dtype=float)
        self.loc = np.array(self.n_veh, dtype=np.array)
        self.load = np.zeros(self.n_veh)
        self.unserved = self.sets['P']   ### Ordered by time window

    def closest_dropoff(self, v):
        closest_drop = None
        min_time = float('inf')  # Initialize with infinity to ensure any valid time will be smaller
        tij = self.params['tij']

        for req in self.vehicle_current_requests[v]:
            drop_node = self.params['pair_pi_di'][req]
            current_time = tij[self.loc[v], drop_node]

            if current_time < min_time:
                min_time = current_time
                closest_drop = drop_node

        return closest_drop

    def request_greedy_heuristic(self):
        tij = self.params['tij']
        while self.unserved:
            p = self.unserved.pop()
            closest = None
            for v in range(self.n_veh):
                if tij[self.loc[v], p] <  tij[self.loc[v], self.vehicle_closest_dropoff[v]]:
                    closest = v
                else:
                    closest = self.vehicle_closest_dropoff[v]
            
            








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


        
    
