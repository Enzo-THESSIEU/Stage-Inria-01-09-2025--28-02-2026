##### Imports #####
import gurobipy as gb
from Parametres_HallPosada import DARPDataBuilder_HallPosada
from Constraints_HallPosada import DARPConstraintBuilder_HallPosada


class DARPModelBuilder_HallPosada:
    """
    Full model builder for the Häll & Posada (2017) Model 2 formulation.
    Homogeneous fleet – single resource.
    """

    def __init__(self, model_name="IDARP_TT_Model2", TIME_LIMIT=600, sets=None, params=None, verbose=True):
        self.model_name = model_name
        self.TIME_LIMIT = TIME_LIMIT
        self.verbose = verbose

        # === Load data ===
        if sets is None or params is None:
            data_builder = DARPDataBuilder_HallPosada()
            self.sets, self.params = data_builder.build()
        else:
            self.sets, self.params = sets, params

        # === Validation ===
        self._validate_params()

        # === Placeholders ===
        self.model = None
        self.vars = {}

    # ----------------------------------------------------------------------
    # Internal validation
    # ----------------------------------------------------------------------
    def _validate_params(self):
        required_keys = ["Cij", "Cij_dr", "Tij", "TD", "TA", "T0j", "Ti_e", "Ti_l", "Bi", "rbar"]
        missing = [k for k in required_keys if k not in self.params]
        if missing:
            raise KeyError(f"❌ Missing required parameters: {missing}")

        for key in ["N", "A", "K", "R"]:
            if key not in self.sets:
                raise KeyError(f"❌ Missing required set: {key}")

    # ----------------------------------------------------------------------
    # Build the MILP model
    # ----------------------------------------------------------------------
    def build(self):
        if self.verbose:
            print(f"🧩 Building model '{self.model_name}'...")

        m = gb.Model(self.model_name)
        m.Params.TimeLimit = self.TIME_LIMIT
        m.Params.OutputFlag = 1 if self.verbose else 0

        # === Constraint Builder ===
        cb = DARPConstraintBuilder_HallPosada(
            m=m,
            vars_={},
            sets=self.sets,
            params=self.params
        )

        # === Variables ===
        vars_ = cb.define_variables()
        cb.vars_ = vars_

        # === Objective ===
        cb.set_objective()

        # === Constraints ===
        cb.add_vehicle_constraints()
        cb.add_passenger_balance_constraints()
        cb.add_capacity_constraints()
        cb.add_time_constraints()
        cb.add_transfer_node_constraints()
        cb.add_strengthening_constraints()

        # === Finalize ===
        cb.finalize()

        self.model, self.vars = m, vars_
        if self.verbose:
            print("✅ Model built successfully.")
        return m, vars_

    # ----------------------------------------------------------------------
    # Run solver
    # ----------------------------------------------------------------------
    def run(self):
        """Run Gurobi optimization and print summary."""
        if self.model is None:
            self.build()

        if self.verbose:
            print("🚀 Optimizing the Hall & Posada (2017) Model 2 formulation...")

        self.model.optimize()

        # === Status handling ===
        status = self.model.status
        if status == gb.GRB.OPTIMAL:
            print(f"✅ Optimal objective value: {self.model.objVal:.2f}")
        elif status == gb.GRB.TIME_LIMIT:
            print(f"⏳ Time limit reached. Best bound: {self.model.ObjBound:.2f}")
        elif status in [gb.GRB.INFEASIBLE, gb.GRB.INF_OR_UNBD]:
            print("❌ Model infeasible or unbounded. Run m.computeIIS() for diagnosis.")
        else:
            print(f"⚠️ Optimization ended with status: {status}")

        return self.model

    # ----------------------------------------------------------------------
    # Extract results summary
    # ----------------------------------------------------------------------
    def summarize_solution(self):
        """Quick textual summary of main variables after optimization."""
        if self.model is None or self.model.status != gb.GRB.OPTIMAL:
            print("⚠️ No optimal solution to summarize.")
            return

        x = self.vars.get("x", {})
        t = self.vars.get("t", {})

        print("\n=== Vehicle Arcs Used (x) ===")
        for key, var in x.items():
            if var.X > 0.5:
                print(f"x{key} = 1")

        print("\n=== Node Service Times (t) ===")
        for node, var in t.items():
            print(f"t{node} = {var.X:.1f}")

# ----------------------------------------------------------------------
# Example Run
# ----------------------------------------------------------------------
if __name__ == "__main__":
    builder = DARPModelBuilder_HallPosada(
        model_name="Hall_Posada_Model2",
        TIME_LIMIT=300,
        verbose=True
    )

    model, vars_ = builder.build()
    builder.run()
    builder.summarize_solution()
