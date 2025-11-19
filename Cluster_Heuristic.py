##### Imports #####
import gurobipy as gb
import random as r

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
    def sym_w_v_inflow(self, i, j):
        """Symmetric connectivity measure for variable-substitution models."""
        A = self.sets['A']
        v = self.vars_['v']
        w = 0.0
        try:
            if (j, i) in v.keys():
                w += self.m.cbGetNodeRel(v[j, i])
        except gb.GurobiError:
            # fallback to current solution value or 0 if not available
            if (j, i) in v.keys():
                w += getattr(v[j, i], "X", 0)
        return w

    def sym_w_v_outflow(self, i, j):
        """Symmetric connectivity measure for variable-substitution models."""
        A = self.sets['A']
        v = self.vars_['v']
        w = 0.0
        try:
            if (i, j) in v.keys():
                w += self.m.cbGetNodeRel(v[i, j])
        except gb.GurobiError:
            # fallback to current solution value or 0 if not available
            if (i, j) in v.keys():
                w += getattr(v[i, j], "X", 0)
        return w


    def sym_w_x_inflow(self, i, j):
        """Symmetric connectivity measure for vehicle-indexed models."""
        A = self.sets['A']
        K = self.sets['K']
        x = self.vars_['x']
        w = 0.0
        try:
            for k in K:
                if (k, j, i) in x.keys():
                    w += self.m.cbGetNodeRel(x[k, j, i])
        except gb.GurobiError:
            # fallback to current or last-known solution
            for k in K:
                if (k, j, i) in x.keys():
                    w += getattr(x[k, j, i], "X", 0)
        return w

    def sym_w_x_outflow(self, i, j):
        """Symmetric connectivity measure for vehicle-indexed models."""
        A = self.sets['A']
        K = self.sets['K']
        x = self.vars_['x']
        w = 0.0
        try:
            for k in K:
                if (k, i, j) in x.keys():
                    w += self.m.cbGetNodeRel(x[k, i, j])
        except gb.GurobiError:
            # fallback to current or last-known solution
            for k in K:
                if (k, i, j) in x.keys():
                    w += getattr(x[k, i, j], "X", 0)
        return w

    # === Clustering heuristic ===
    def cluster_builder(self, max_weight=1.0):
        """
        Builds clusters (connected components) based on current integer solution.
        """
        if self.variable_substitution:
            sym_func_inflow  =  self.sym_w_v_inflow
            sym_func_outflow = self.sym_w_v_outflow
        else:
            sym_func_inflow  =  self.sym_w_x_inflow
            sym_func_outflow = self.sym_w_x_outflow

        P_D = self.sets['P'] + self.sets['D']
        N = self.sets['N']
        A = self.sets['A']
        zeroDepot_node, endDepot_node = self.sets['zeroDepot'],  self.sets['endDepot']

        clusters = []

        # Only consider non-depot nodes
        candidate_nodes = {n for n in N if n not in [zeroDepot_node, endDepot_node]}

        for i in P_D:

            cluster = {i}
            rest_nodes = candidate_nodes - cluster

            inflow = 0.0
            outflow = 0.0

            expanded = True
            while expanded:

                if rest_nodes:

                    # pick candidates that have at least one connection to the cluster
                    connected_candidates = [
                        j for j in rest_nodes
                        if any((j, c) or (c, j) in A for c in cluster)
                    ]

                    if connected_candidates:
                        random_node = r.choice(connected_candidates)
                    else:
                        random_node = None

                    cluster.add(random_node)

                    rest_nodes = candidate_nodes - cluster
                    inflow = sum(sym_func_inflow(rest_node , cluster_node) for rest_node in rest_nodes for cluster_node in cluster)
                    outflow = sum(sym_func_outflow(cluster_node,rest_node) for rest_node in rest_nodes for cluster_node in cluster)

                    if inflow < max_weight or outflow < max_weight:
                        clusters.append(cluster)
                        expanded = False
                else:
                    expanded = False

        return clusters
