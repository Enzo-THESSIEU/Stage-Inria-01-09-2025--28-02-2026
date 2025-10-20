##### Imports #####
import gurobipy as gb

class DARPHeuristic:
    """
    Callback-based heuristic and subtour elimination cuts for DARP/IDARP models.
    Can be attached to Gurobi via m.optimize(self.subtour_callback).
    """

    def __init__(self, m, vars_, sets, params, variable_substitution=True):
        self.variable_substitution = variable_substitution
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params

    # === Helper weight access ===
    def sym_w_v(self, i, j):
        """Symmetric connectivity measure for variable-substitution models."""
        A = self.sets['A']
        v = self.vars_['v']
        w = 0.0
        try:
            if (i, j) in A:
                w += self.m.cbGetNodeRel(v[i, j])
            if (j, i) in A:
                w += self.m.cbGetNodeRel(v[j, i])
        except gb.GurobiError:
            # fallback to current solution value or 0 if not available
            if (i, j) in A:
                w += getattr(v[i, j], "X", 0)
            if (j, i) in A:
                w += getattr(v[j, i], "X", 0)
        return w


    def sym_w_x(self, i, j):
        """Symmetric connectivity measure for vehicle-indexed models."""
        A = self.sets['A']
        K = self.sets['K']
        x = self.vars_['x']
        w = 0.0
        try:
            if (i, j) in A:
                for k in K:
                    w += self.m.cbGetNodeRel(x[k, i, j])
            if (j, i) in A:
                for k in K:
                    w += self.m.cbGetNodeRel(x[k, j, i])
        except gb.GurobiError:
            # fallback to current or last-known solution
            if (i, j) in A:
                for k in K:
                    w += getattr(x[k, i, j], "X", 0)
            if (j, i) in A:
                for k in K:
                    w += getattr(x[k, j, i], "X", 0)
        return w


    # === Clustering heuristic ===
    def cluster_builder(self, max_weight=1.0):
        """
        Builds clusters (connected components) based on current integer solution.
        """
        if self.variable_substitution:
            sym_func = self.sym_w_v
        else:
            sym_func = self.sym_w_x

        P_D = self.sets['P'] + self.sets['D']
        N = self.sets['N']
        zeroDepot_node, endDepot_node = self.sets['zeroDepot'],  self.sets['endDepot']

        clusters = []
        visited = set()

        # Only consider non-depot nodes
        candidate_nodes = {n for n in N if n not in [zeroDepot_node, endDepot_node]}

        for i in P_D:
            if i in visited:
                continue  # skip nodes already assigned to some cluster

            cluster = {i}
            rest_nodes = candidate_nodes - cluster

            inflow = 0.0
            outflow = 0.0

            expanded = True
            while expanded:
                expanded = False
                new_nodes = set()

                # Check all remaining nodes for potential connection
                for j in rest_nodes:
                    # Calculate connectivity contribution with current cluster
                    inflow_j = sum(sym_func(j, c) for c in cluster)
                    outflow_j = sum(sym_func(c, j) for c in cluster)
                    inflow += inflow_j
                    outflow += outflow_j

                    # Condition to include j: inflow/outflow still below threshold
                    if inflow + outflow < max_weight:
                        new_nodes.add(j)

                if new_nodes:
                    cluster |= new_nodes
                    rest_nodes -= new_nodes
                    expanded = True

                # Check if we’ve reached a strong enough cluster
                if inflow + outflow >= max_weight:
                    break

            # Store this cluster
            clusters.append(cluster)
            visited |= cluster
        # print("clusters:",clusters)

        return clusters
