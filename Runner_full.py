##### Imports #####
from tqdm import tqdm
import time as t
import itertools
import gurobipy as gb
import pandas as pd

from Model import DARPModelBuilder
from routes_2 import extract_route_final
from Parametres import DARPDataBuilder
from Cluster_Heuristic import DARPHeuristic  # assuming you renamed Cluster_Heuristic to this

class DARPExperimentRunner:
    def __init__(self, time_limit=2 * 60 * 60):
        """Initialize experiment runner."""
        self.time_limit = time_limit
        self.results = []
        self.skipped = []

    def stop_when_optimal(self, model, where):
        """Stop solver when objective 880 is found, and record time."""
        if where == gb.GRB.Callback.MIPSOL:
            obj_val = model.cbGet(gb.GRB.Callback.MIPSOL_OBJ)
            if abs(obj_val - 880) < 1e-6:  # tolerance for numerical precision
                model._t_optimal = t.perf_counter() - model._start_time
                print(f"\n🎯 Optimal solution (880) found after {model._t_optimal:.2f} seconds — stopping solver.")
                model.terminate()

    # === SINGLE MODEL RUN ===
    def run_single_model(self, bool_params):
        """Run a single configuration of the DARP model."""
        model_name = "IDARP_Model_Single"
        print(f"\n=== Running single model: {model_name} ===")
        print(bool_params)

        # Build data
        data_builder = DARPDataBuilder(
            duplicate_transfers=bool_params["duplicate_transfers"],
            arc_elimination=bool_params["arc_elimination"],
            ev_constraints=bool_params["ev_constraints"],
            use_imjn=bool_params["use_imjn"]
        )
        data_sets, data_params = data_builder.build()

        clusters = []

        # Build model
        model_builder = DARPModelBuilder(
            model_name=model_name,
            TIME_LIMIT=self.time_limit,
            sets=data_sets,
            params=data_params,
            options = bool_params
        )

        m, vars_ = model_builder.build()

        # Attach callback data
        m._v = vars_["v"]
        m._x = vars_["x"]
        m._P = data_sets["P"]
        m._D = data_sets["D"]
        m._N = data_sets["N"]
        m._K = data_sets["K"]
        m._A = data_sets["A"]
        m.Params.LazyConstraints = 1

        # === Optimization ===
        start_cpu = t.process_time()
        start_wall = t.perf_counter()

        m._start_time = start_wall      # start reference for callback
        m._t_optimal = None             # initialize storage variable

        # --- Run with subtour elimination or not ---
        if bool_params['subtour_elimination']:
            heuristic = DARPHeuristic(
                m, vars_, data_sets, data_params,
                variable_substitution=bool_params['variable_substitution']
            )

            def lp_callback(model, where):
                if where == gb.GRB.Callback.MIPNODE:
                    clusters = heuristic.cluster_builder(max_weight=1.0)

                    # access sets and variables
                    A = heuristic.sets["A"]
                    N = heuristic.sets["N"]
                    K = heuristic.sets["K"]
                    v = heuristic.vars_["v"]
                    x = heuristic.vars_["x"]

                    for cluster in clusters:
                        if len(cluster) <= 1:
                            continue  # skip singletons

                        if heuristic.variable_substitution:
                            expr = gb.quicksum(
                                v[i, j]
                                for i in cluster
                                for j in N
                                if j not in cluster and (i, j) in A
                            )
                        else:
                            expr = gb.quicksum(
                                x[k, i, j]
                                for k in K
                                for i in cluster
                                for j in N
                                if j not in cluster and (i, j) in A
                            )

                        # === ADD CUT to current LP relaxation ===
                        model.cbCut(expr >= 1)
                        print(f"Added user cut for cluster {cluster}")
                self.stop_when_optimal(model, where)

            m.optimize(lp_callback)
        else:
            m.optimize(lambda model, where: self.stop_when_optimal(model, where))

        end_cpu = t.process_time()
        end_wall = t.perf_counter()
        cpu_time = end_cpu - start_cpu
        elapsed_time = end_wall - start_wall

        print("\n=== Vehicle arc variables (x[k,i,j]) with nonzero values ===")
        for (k, i, j), var in vars_['x'].items():
            try:
                val = var.X
            except gb.GurobiError:
                val = float('nan')
            if abs(val) > 1e-6:
                print(f"x[{k},{i},{j}] = {val:.3f}")



        # === Handle infeasibility ===
        if m.status == gb.GRB.INFEASIBLE:
            print("⚠️ Model infeasible — writing IIS files...")
            m.computeIIS()
            m.write(f"{model_name}.ilp")
            m.write(f"{model_name}.iis")
            raise SystemExit("Model infeasible; IIS written.")

        # === Extract routes ===
        v1, v2, r1, r2, r3, r4 = extract_route_final(
            vars_,
            nodes=data_sets["nodes"],
            N=data_sets["N"],
            K=data_sets["K"],
            variable_substitution=bool_params['variable_substitution']
        )

        # === Save results ===
        result = {
            "Model": model_name,
            **bool_params,
            "Objective": m.ObjVal,
            "Lower Bound": m.ObjBound,
            "Vehicle 1 Route": str(v1),
            "Vehicle 2 Route": str(v2),
            "Request 1 Route": str(r1),
            "Request 2 Route": str(r2),
            "Request 3 Route": str(r3),
            "Request 4 Route": str(r4),
            "CPU Time (s)": cpu_time,
            "Elapsed Time (s)": elapsed_time,
            
        }
        self.results.append(result)
        print(result)
        print(f"\n✅ Done: {model_name} | Obj = {m.ObjVal:.2f}, Bound = {m.ObjBound:.2f}")
        return result

    # === MULTIPLE MODEL RUN ===
    def run_all_combinations(self):
        """Run all combinations of Boolean model parameters."""
        boolean_params = {
            'roaund': [True, False],
            'robund': [True, False],
            # 'duplicate_transfers': [True, False],
            'arc_elimination': [True, False],
            'variable_substitution': [True, False],
            'subtour_elimination': [True, False],
            # 'timetabled_departures': [True, False],
            # 'use_imjn': [True, False],
        }                 

        all_combos = list(itertools.product(*boolean_params.values()))
        total = len(all_combos)
        print(f"=== Starting {total} model configurations ===")

        pbar = tqdm(total=total, desc="Progress", unit="config", ncols=100)
        start_global = t.perf_counter()

        for idx, combo in enumerate(all_combos, start=1):
            params = dict(zip(boolean_params.keys(), combo))
            model_name = "IDARP_Model_" + "_".join([f"{k[:3]}{int(v)}" for k, v in params.items()])

            # Default modifiers
            params['ev_constraints'] = False
            params['transfer_node_strengthening'] = True
            params['timetabled_departures']= False
            params['use_imjn']= False
            params['duplicate_transfers']= True

            # Feasibility filters
            skip_reason = self.check_invalid_combo(params)
            if skip_reason:
                self.skipped.append({**params, "Reason": skip_reason})
                pbar.update(1)
                pbar.set_postfix({"Status": "⏭️ skipped"})
                continue

            # === Build model ===
            data_builder = DARPDataBuilder(
                duplicate_transfers=params["duplicate_transfers"],
                arc_elimination=params["arc_elimination"],
                ev_constraints=params["ev_constraints"],
                use_imjn=params["use_imjn"]
            )
            data_sets, data_params = data_builder.build()


            model_builder = DARPModelBuilder(
                model_name=model_name,
                TIME_LIMIT=self.time_limit,
                sets=data_sets,
                params=data_params,
                options=params
            )
            m, vars_ = model_builder.build()

            heuristic = DARPHeuristic(m, vars_, data_sets, data_params,
                                      variable_substitution=params["variable_substitution"])

            # Attach callback data
            m._v = vars_["v"]
            m._x = vars_["x"]
            m._P = data_sets["P"]
            m._D = data_sets["D"]
            m._N = data_sets["N"]
            m._K = data_sets["K"]
            m._A = data_sets["A"]
            m.Params.LazyConstraints = 1

            # # === Solve ===
            # start_wall = t.perf_counter()
            # if params['subtour_elimination']:
            #     m.optimize(heuristic.subtour_callback)
            # else:
            #     m.optimize()
            # end_wall = t.perf_counter()
            # elapsed_time = end_wall - start_wall

            # === Optimization ===
            start_cpu = t.process_time()
            start_wall = t.perf_counter()

            m._start_time = start_wall      # start reference for callback
            m._t_optimal = None             # initialize storage variable

            # --- Run with subtour elimination or not ---
            if params['subtour_elimination']:
                heuristic = DARPHeuristic(
                    m, vars_, data_sets, data_params,
                    variable_substitution=params['variable_substitution']
                )

                def combined_callback(model, where):
                    # Apply subtour cuts and stop when optimal found
                    heuristic.subtour_callback(model, where)
                    self.stop_when_optimal(model, where)

                m.optimize(combined_callback)
            else:
                m.optimize(lambda model, where: self.stop_when_optimal(model, where))

            end_cpu = t.process_time()
            end_wall = t.perf_counter()
            cpu_time = end_cpu - start_cpu
            elapsed_time = end_wall - start_wall

            # === Record ===
            result = self._record_result(m, params, data_sets, vars_, elapsed_time)
            self.results.append(result)

            # Update progress
            avg_time = (t.perf_counter() - start_global) / idx
            remaining = avg_time * (total - idx)
            eta_min, eta_sec = divmod(int(remaining), 60)
            pbar.update(1)
            pbar.set_postfix({"ETA": f"{eta_min:02d}:{eta_sec:02d}", "Last": f"{elapsed_time:.1f}s"})

        pbar.close()
        self._save_results()

    # === UTILITIES ===
    def check_invalid_combo(self, params):
        """Return reason string if combo invalid, else None."""
        if params["timetabled_departures"] and not params["use_imjn"]:
            return "Timetabled departures require IMJN nodes (use_imjn=True)."
        if params["use_imjn"] and params["duplicate_transfers"]:
            return "Duplicate transfers cannot be used with IMJN nodes."
        if params["timetabled_departures"] and not params["variable_substitution"]:
            return "Timetabled Departures without variable substitution are inconsistent."
        if not params["duplicate_transfers"] and not params["use_imjn"]:
            return "Cannot disable duplicate_transfers without IMJN nodes."
        return None

    def _record_result(self, m, params, data_sets, vars_, elapsed_time):
        """Record one experiment’s result row."""
        if m.status == gb.GRB.INFEASIBLE:
            print(f"❌ {params} infeasible — skipping.")
            return {
                "Model": "Infeasible",
                **params,
                "Objective": "Infeasible",
                "Lower Bound": "N/A",
                "Elapsed Time (s)": elapsed_time,
            }
        elif m.SolCount == 0:
            print(f"⚠️ {params} — No feasible solution or timeout.")
            return {
                "Model": "IDARP_Model_" + "_".join(
                    [f"{k[:3]}{int(v)}" for k, v in params.items() 
                    if k in ["roaund", "robund", "arc_elimination", "variable_substitution", 
                            "subtour_elimination"]]
                ),
                **params,
                "Objective": "No feasible solution",
                "Lower Bound": "N/A",
                "Elapsed Time (s)": elapsed_time,
            }
        else:
            print(f"✅ Solved — Obj = {m.ObjVal:.2f}, Bound = {m.ObjBound:.2f}")
            try:
                v1, v2, r1, r2, r3, r4 = extract_route_final(
                    vars_,
                    nodes=data_sets['nodes'],
                    N=data_sets['N'],
                    K=data_sets['K'],
                    variable_substitution=params["variable_substitution"]
                )
            except Exception as e:
                print(f"⚠️ Route extraction failed: {e}")
                v1 = v2 = r1 = r2 = r3 = r4 = []

            return {
                "Model" : "IDARP_Model_" + "_".join(
                    [f"{k[:3]}{int(v)}" for k, v in params.items() 
                    if k in ["roaund", "robund", "arc_elimination", "variable_substitution", 
                            "subtour_elimination"]]
                ),
                **params,
                "Objective": m.ObjVal,
                "Lower Bound": m.ObjBound,
                "Vehicle 1 Route": str(v1),
                "Vehicle 2 Route": str(v2),
                "Request 1 Route": str(r1),
                "Request 2 Route": str(r2),
                "Request 3 Route": str(r3),
                "Request 4 Route": str(r4),
                "Elapsed Time (s)": elapsed_time,
            }

    def _save_results(self):
        """Save results and skipped combinations to CSV."""
        df = pd.DataFrame(self.results)
        df.to_csv("results_full.csv", index=False)
        print("\n✅ All feasible results saved to results_full.csv")

        if self.skipped:
            df_skipped = pd.DataFrame(self.skipped)
            df_skipped.to_csv("skipped_combinations.csv", index=False)
            print(f"⚙️ Skipped {len(self.skipped)} invalid configurations → saved to skipped_combinations.csv")


# === Example usage ===
if __name__ == "__main__":
    TIME_LIMIT = 2 * 60 * 60
    bool_params_singular = {
        "duplicate_transfers": True,
        "arc_elimination": True,
        "variable_substitution": True,
        "subtour_elimination": True,
        "transfer_node_strengthening": True,
        "ev_constraints": False,
        "timetabled_departures": False,
        "use_imjn": False,
    }

    runner = DARPExperimentRunner(time_limit=TIME_LIMIT)
    runner.run_single_model(bool_params_singular)
    # runner.run_all_combinations()
