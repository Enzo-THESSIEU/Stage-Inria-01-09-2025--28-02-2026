class DARPDebuggingFunctions:
    """Utility class for debugging and analyzing DARP model outputs."""

    def __init__(self, m, vars_, sets, params):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params
        self.base = lambda i: i[0] if isinstance(i, tuple) else i

    def extract_variable_values(self):
        """Extract and return all main variable values."""
        z_values_all = {k: int(v.X) for k, v in self.vars_['z'].items()}
        z_values_1 = {k: int(v.X) for k, v in self.vars_['z'].items() if abs(v.X) > 1e-6}

        y_values_r4 = {k: int(v.X) for k, v in self.vars_['y'].items() if k[0] == 4}
        a_values = {k: int(v.X) for k, v in self.vars_['a'].items()}

        v_values_1 = {
            k: self.params['tij'][k]
            for k, v in self.vars_['v'].items()
            if abs(v.X) > 1e-6
        }
        sum_v_values = sum(v_values_1.values())

        return z_values_all, z_values_1, y_values_r4, a_values, v_values_1, sum_v_values

    def check_transfer_balance(self, extractor):
        """Check PT transfer balance using extractor."""
        sum_balance = extractor.test_passenger_transfer_balance(*extractor.extract_vehicle_route_final())
        wrong_sum_balance = {
            k: v for k, v in sum_balance.items() if abs(v['diff']) > 1e-6
        }
        return wrong_sum_balance

    def compute_constraint_balance(self):
        """Compute passenger flow balance for each request and node."""
        R, C, N, A = self.sets['R'], self.sets['C'], self.sets['N'], self.sets['A']
        Departures = self.params['Departures']
        fi_r = self.params['fi_r']

        Constraint_val = {}
        for r in R:
            for i in C:
                # Passenger movements
                sum_pickup = sum(
                    self.vars_['y'][r, i, j].X for j in N if (i, j) in A and (r, i, j) in self.vars_['y']
                )
                sum_dropoff = sum(
                    self.vars_['y'][r, j, i].X for j in N if (j, i) in A and (r, j, i) in self.vars_['y']
                )

                # PT transfers
                z_sum_pickup, z_sum_dropoff = 0, 0
                for j in C:
                    if (self.base(i), self.base(j)) in Departures:
                        for d in Departures[(self.base(i), self.base(j))]:
                            if (d, i, j) in self.vars_['z']:
                                z_sum_pickup += self.vars_['z'][(d, i, j)].X
                    if (self.base(j), self.base(i)) in Departures:
                        for d in Departures[(self.base(j), self.base(i))]:
                            if (d, j, i) in self.vars_['z']:
                                z_sum_dropoff += self.vars_['z'][(d, j, i)].X

                lhs = sum_pickup + z_sum_pickup - sum_dropoff - z_sum_dropoff
                rhs = fi_r.get((r, i), 0)

                Constraint_val[(r, i)] = {
                    'y_out': sum_pickup,
                    'y_in': sum_dropoff,
                    'z_out': z_sum_pickup,
                    'z_in': z_sum_dropoff,
                    'lhs': lhs,
                    'rhs': rhs,
                    'diff': lhs - rhs,
                }
        return Constraint_val
    
    def print_z_variable_summary(self, model_name):
        """
        Print a summary of z variables, check for missing PT arcs,
        and print final objective info.
        """
        vars_ = self.vars_
        C = self.sets["C"]
        Departures = self.params["Departures"]

        Total_z_variables = len(vars_['z'])
        first_ten_elements = list(vars_['z'].keys())[:10]

        count_existing = 0
        count_missing = 0

        for i in C:
            for j in C:
                if (self.base(i), self.base(j)) in Departures:
                    for d in Departures[(self.base(i), self.base(j))]:
                        if (d, i, j) in vars_['z']:
                            count_existing += 1
                        else:
                            count_missing += 1

        return Total_z_variables, first_ten_elements, count_existing, count_missing
