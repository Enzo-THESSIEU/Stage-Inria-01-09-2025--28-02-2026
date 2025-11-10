import gurobipy as gb
from Parametres import DARPDataBuilder
from Constraints import DARPConstraintBuilder

### Rewritten as classes

class DARPModelBuilder:
    def __init__(self, model_name, options, TIME_LIMIT, sets=None, params=None):
        self.model_name = model_name
        self.options = options
        self.TIME_LIMIT = TIME_LIMIT

        if sets is None or params is None:
            data_builder = DARPDataBuilder(
                duplicate_transfers=options.get("duplicate_transfers", True),
                arc_elimination=options.get("arc_elimination", True),
                ev_constraints=options.get("ev_constraints", False),
                use_imjn=options.get("use_imjn", False)
            )
            self.sets, self.params = data_builder.build()
        else:
            self.sets, self.params = sets, params

        self.model = None
        self.vars = {}

    def build(self):
        # Unpack options
        duplicate_transfers = self.options.get("duplicate_transfers", True)
        arc_elimination = self.options.get("arc_elimination", True)
        variable_substitution = self.options.get("variable_substitution", True)
        subtour_elimination = self.options.get("subtour_elimination", True)
        transfer_node_strengthening = self.options.get("transfer_node_strengthening", True)
        ev_constraints = self.options.get("ev_constraints", False)
        timetabled_departures = self.options.get("timetabled_departures", False)
        use_imjn = self.options.get("use_imjn", False)
        MoPS = self.options.get("MoPS", False)
        clusters = []

        # === Step 2: Build Model ===
        m = gb.Model(self.model_name)
        m.setParam("TimeLimit", self.TIME_LIMIT)

        # === Step 3: Create Constraint Builder ===
        constraint_builder = DARPConstraintBuilder(
            m=m,
            vars_={},  # will be populated below
            sets=self.sets,
            params=self.params,
            duplicate_transfers = duplicate_transfers,
            arc_elimination = arc_elimination,
            variable_substitution = variable_substitution,
            subtour_elimination = subtour_elimination,
            transfer_node_strengthening = transfer_node_strengthening,
            ev_constraints = ev_constraints,
            timetabled_departures = timetabled_departures,
            use_imjn = use_imjn
        )

        # === Step 4: Define Variables ===
        vars_ = constraint_builder.define_variables()
        constraint_builder.vars_ = vars_  # store them back for other methods

        # === Step 5: Objective ===
        constraint_builder.variable_substitution = variable_substitution
        constraint_builder.timetabled_departures = timetabled_departures
        constraint_builder.set_objective()

        # === Step 6: Constraints ===
        constraint_builder.add_vehicle_logic_constraints()
        constraint_builder.add_passenger_balance_constraints()
        constraint_builder.add_capacity_constraints()
        constraint_builder.add_time_consistency_constraints()

        if variable_substitution:
            constraint_builder.add_variable_substitution_constraints()

        if transfer_node_strengthening:
            constraint_builder.add_transfer_node_constraints()

        if ev_constraints:
            constraint_builder.add_battery_constraints()

        if timetabled_departures:
            constraint_builder.add_scheduled_PT_constraints()

        if use_imjn:
            op = 1
            constraint_builder.add_artificial_node_constraints()

        if MoPS:
            constraint_builder.add_MoPS_constraints()

        # constraint_builder.base_model_optimal_solution_constrainer_paper()

        constraint_builder.debugging_constraint()

        # === Finalize ===
        m.update()
        self.model, self.vars = m, vars_
        return self.model, self.vars
