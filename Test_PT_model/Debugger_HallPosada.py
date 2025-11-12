class DARPFlowDebugger_HallPosada:
    """
    Debug utility to compute passenger flow balance at every node for each request.
    Works with solved Gurobi model and variable dictionaries (y, z).
    """

    def __init__(self, model, vars_, sets, params):
        self.m = model
        self.vars_ = vars_
        self.sets = sets
        self.params = params

    def base(self, x):
        """Return flat node id from (i,m) or int."""
        return x[0] if isinstance(x, tuple) else x

    def check_passenger_flow_balance(self, tol=1e-5):
        """
        Computes inflow, outflow and balance of passenger flows for each node i and each request r.
        Prints imbalances larger than tolerance.
        """
        y = self.vars_.get("y", {})
        z = self.vars_.get("z", {})
        A = self.sets.get("A", [])
        R = self.sets.get("R", [])
        NG = self.sets.get("NG", [])
        NP = self.sets.get("NP", [])
        ND = self.sets.get("ND", [])
        Fir = self.params.get("Fir", {})

        print("\n=== DEBUG: Passenger Flow Balance by Node and Request ===")
        for r in R:
            print(f"\nRequest {r}:")
            for node in NP + ND + NG:
                i = node
                inflow_y = sum(y[r, j, i].X for (j, jj) in A if jj == i and (r, j, i) in y)
                outflow_y = sum(y[r, i, j].X for (ii, j) in A if ii == i and (r, i, j) in y)

                # Include PT inflows/outflows if PT arcs exist
                inflow_z = sum(z[jj, self.base(i), d, r].X
                               for (jj, ii, d, rr) in z.keys()
                               if ii == self.base(i) and rr == r)
                outflow_z = sum(z[self.base(i), jj, d, r].X
                                for (ii, jj, d, rr) in z.keys()
                                if ii == self.base(i) and rr == r)

                rhs = Fir.get((r, i), 0)
                balance = (outflow_y + outflow_z) - (inflow_y + inflow_z)
                gap = balance - rhs

                if abs(gap) > tol:
                    print(f"⚠️ Node {i}: inflow={inflow_y + inflow_z:.2f}, "
                          f"outflow={outflow_y + outflow_z:.2f}, "
                          f"rhs={rhs}, balance={balance:.2f}, gap={gap:.2e}")
                else:
                    print(f"✓ Node {i}: inflow={inflow_y + inflow_z:.2f}, "
                          f"outflow={outflow_y + outflow_z:.2f}, "
                          f"rhs={rhs}, OK")

        print("\n=== End of Passenger Flow Balance Debug ===")
