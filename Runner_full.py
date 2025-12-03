##### Imports #####
from tqdm import tqdm
import time as t
import itertools
import gurobipy as gb
import pandas as pd
import ast

from Model import DARPModelBuilder
from routes_2 import DARPRouteExtractor
from Parametres import DARPDataBuilder
from Cluster_Heuristic import DARPHeuristic 
from Debugger_code import DARPDebuggingFunctions, DARPRouteDebugging
from model_verification import write_model_verification_report
from LNS import initial_solution, translate_LNS_to_Gurobi

import logging
from logging.handlers import RotatingFileHandler
import json
from pathlib import Path
from datetime import datetime

class DARPLogger:
    def __init__(self, log_dir="logs", json_enabled=True):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(exist_ok=True)

        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self.log_file = self.log_dir / f"darp_run_{timestamp}.log"
        self.json_file = self.log_dir / f"darp_run_{timestamp}.json"

        # --- Setup text logger ---
        self.logger = logging.getLogger("DARPLogger")
        self.logger.setLevel(logging.INFO)

        handler = RotatingFileHandler(self.log_file, maxBytes=5_000_000, backupCount=5)
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S"
        )
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

        self.json_enabled = json_enabled
        if self.json_enabled:
            self.json_data = []

    def log(self, message):
        """Write simple text log."""
        self.logger.info(message)

    def log_dict(self, data_dict, title=None):
        """Text log + JSON structured record for results."""
        if title:
            self.logger.info(f"=== {title} ===")
        self.logger.info(json.dumps(data_dict, indent=2))

        if self.json_enabled:
            entry = {"title": title, "data": data_dict}
            self.json_data.append(entry)

    def save_json(self):
        """Write accumulated JSON logs."""
        if self.json_enabled:
            with open(self.json_file, "w") as f:
                json.dump(self.json_data, f, indent=2)



class DARPExperimentRunner:

    def __init__(self, time_limit=2*60*60, stream_file="results_stream.csv"):
        self.time_limit = time_limit
        self.stream_file = stream_file

        self.results = []
        self.skipped = []
        self.debugging_data_all = []

        self.logger = DARPLogger()

        # streaming: remove old file if exists
        if Path(self.stream_file).exists():
            Path(self.stream_file).unlink()

    # ---------------------------
    # BASIC HELPERS
    # ---------------------------

    def _append_csv_row(self, row):
        """Write exactly one row incrementally to CSV."""
        df = pd.DataFrame([row])
        header = not Path(self.stream_file).exists()
        df.to_csv(self.stream_file, mode='a', header=header, index=False)

    def _model_id(self, params):
        return "IDARP_" + "_".join([f"{k[:3]}{int(v)}" for k, v in params.items()])

    # ---------------------------
    # 1 — BUILD
    # ---------------------------

    def build_model(self, params, model_name):

        self.logger.log(f"Building model {model_name}")

        # === 1. Build model data (sets + parameters) ===
        data_builder = DARPDataBuilder(
            duplicate_transfers=params["duplicate_transfers"],
            arc_elimination=params["arc_elimination"],
            ev_constraints=params["ev_constraints"],
            use_imjn=params["use_imjn"],
            MoPS=params["MoPS"]
        )

        sets, p = data_builder.build()

        # === 2. Build Gurobi model with all variables & constraints ===
        model_builder = DARPModelBuilder(
            model_name=model_name,
            TIME_LIMIT=self.time_limit,
            sets=sets,
            params=p,
            options=params
        )

        m, vars_ = model_builder.build()

        # Store
        self.m, self.vars_, self.sets, self.params = m, vars_, sets, p

        self.m.update()


        # ==========================================================================
        # === 3. INITIAL SOLUTION → Construct heuristic solution ===================
        # ==========================================================================

        if params['LNS']:

            # Build heuristic solution
            init = initial_solution(
                m=self.m,
                vars_=self.vars_,
                sets=self.sets,
                params=self.params,
                **params
            )

            # Run heuristic of your choice
            # init.closest_request_greedy_heuristic()
            # init.earliest_request_greedy_heuristic()
            init.Execute_Heuristic()
            vehicle_routes = init.vehicle_route

            # === LOG heuristic routes ===
            for k in range(len(vehicle_routes)):
                print(f"Initial route for vehicle {k}: {vehicle_routes[k]}")
                self.logger.log(f"Initial route for vehicle {k}: {vehicle_routes[k]}")


            # ==========================================================================
            # === 4. Translate heuristic into MIP start ================================
            # ==========================================================================

            

            translator = translate_LNS_to_Gurobi(
                vehicle_routes=vehicle_routes,
                m=self.m,
                vars_=self.vars_,
                sets=self.sets,
                params=self.params,
                variable_substitution=params["variable_substitution"],
                duplicate_transfers=params["duplicate_transfers"],
                use_imjn=params["use_imjn"]
            )

            # Set initial x/v values
            translator.intialise_all_gurobi_binary_variables_to_one()
            translator.translate_to_gurobi_variables_x()

            # Set y/z depending on duplicate vs IMJN
            if params["duplicate_transfers"]:
                translator.translate_to_gurobi_variables_y_z_duplicate_transfers()
                print("runnning translator.translate_to_gurobi_variables_y_z_duplicate_transfers()")
            else:
                translator.translate_to_gurobi_variables_y_z_imjn()
                print("running translator.translate_to_gurobi_variables_y_z_imjn()")

            # Push warm-start info into Gurobi
            self.m.update()

            print("\nMIP start objective (if any):", self.m.NumStart)

            print("\n=== CHECKING ASSIGNED START VALUES ===")

            count = 0

            for var in self.m.getVars():
                # Skip z-variables
                if var.VarName.startswith("z"):
                    continue

                # Check if Start value was assigned
                if var.Start is not None and 0.5 < abs(var.Start) < 1e9:
                    print(f"{var.VarName} = {var.Start}")
                    count += 1

            print(f"Total warm-start variables assigned: {count}")



            self.logger.log("Initial solution applied as MIP start")

        return self.m, self.vars_

    # ---------------------------
    # 2 — Test initial solution
    # ---------------------------
    def test_mip_start_feasibility(self, m):
        # Make a copy so you don't touch the real model
        mtest = m.copy()

        # Fix all integer/binary vars to their Start value
        for v in mtest.getVars():
            if v.VType != gb.GRB.CONTINUOUS and v.Start not in (None, 0):
                v.LB = v.Start
                v.UB = v.Start

        mtest.Params.OutputFlag = 1
        mtest.Params.TimeLimit = 60  # small time limit is enough
        print("\n=== Testing warm-start feasibility ===")
        mtest.optimize()

        print("Test status:", mtest.Status)
        if mtest.Status == gb.GRB.INFEASIBLE:
            print("Warm start is infeasible → computing IIS...")
            mtest.computeIIS()
            mtest.write("warmstart_iis.ilp")
            print("IIS written to warmstart_iis.ilp")

    # ---------------------------
    # 3 — OPTIMISE
    # ---------------------------

    def optimise_model(self, params, model_name):

        m = self.m
        m._v = self.vars_["v"]
        m._x = self.vars_["x"]
        m._P = self.sets["P"]
        m._D = self.sets["D"]
        m._N = self.sets["N"]
        m._K = self.sets["K"]
        m._A = self.sets["A"]

        m.Params.LazyConstraints = 1

        start_cpu = t.process_time()
        start_wall = t.perf_counter()

        self.m.Params.StartNumber = 0
        self.m.Params.StartNodeLimit = 1


        if params["subtour_elimination"]:
            heuristic = DARPHeuristic(
                m, self.vars_, self.sets, self.params,
                variable_substitution=params["variable_substitution"]
            )

            m.Params.Presolve = 0
            m.Params.Aggregate = 0
            m.Params.PreSparsify = 0

            def subtour_callback(model, where):
                if where == gb.GRB.Callback.MIPNODE:
                    try:
                        if params["use_imjn"]:
                            _ = model.cbGetNodeRel(model._x[0,(0,1),(9,1)])
                        else:
                            _ = model.cbGetNodeRel(model._x[0,0,9])
                    except gb.GurobiError:
                        return

                    clusters = heuristic.cluster_builder(max_weight=1.0)

                    v = heuristic.vars_["v"]
                    x = heuristic.vars_["x"]
                    N = heuristic.sets["N"]
                    K = heuristic.sets["K"]

                    for C in clusters:
                        if len(C) <= 1:
                            continue

                        if heuristic.variable_substitution:
                            expr = gb.quicksum(
                                v[i,j]
                                for i in C for j in N
                                if j not in C and (i,j) in v
                            )
                        else:
                            expr = gb.quicksum(
                                x[k,i,j]
                                for k in K for i in C for j in N
                                if j not in C and (k,i,j) in x
                            )

                        model.cbCut(expr >= 1)

            # print("MIP start objective:", self.m.MIPStart)

            m.optimize(subtour_callback)

        else:
            # print("MIP start objective:", self.m.MIPStart)

            m.optimize()




        self.cpu_time = t.process_time() - start_cpu
        self.elapsed_time = t.perf_counter() - start_wall

        if m.status == gb.GRB.INFEASIBLE:
            self.logger.log(f"Model {model_name} is infeasible.")
            m.computeIIS()
            m.write("Conflicting_Constraints.ilp")
            raise SystemExit("Model infeasible; IIS written.")

    # ---------------------------
    # 4 — DEBUG / EXTRACT ROUTES
    # ---------------------------

    def debug_model(self, params, model_name):

        debugger = DARPDebuggingFunctions(
            self.m, self.vars_, self.sets, self.params, options=params
        )

        # constrained, problems = debugger.compute_constraint_balance()
        # vflow, vprob = debugger.check_flow_conservation()
        # debugger.print_z_variable_summary(model_name)

        routedbg = DARPRouteDebugging(
            m=self.m, vars_=self.vars_, sets=self.sets, params=self.params, options=params
        )

        if params['ev_constraints']:
            v1, v2 = routedbg.extract_vehicle_route_final_ev()
            ev_constraint = routedbg.ev_constraints_issue()
        else:
            v1, v2 = routedbg.extract_vehicle_route_final()
            ev_constraint = None
        
        r1, r2, r3, r4 = routedbg.extract_request_route_final()
        z1 = routedbg.extract_PT_route_final()

        self.result = {
            "Model": model_name,
            **params,
            "Objective": self.m.ObjVal,
            "Lower Bound": self.m.ObjBound,
            "CPU Time": self.cpu_time,
            "Wall Time": self.elapsed_time,
            "Veh1": str(v1),
            "Veh2": str(v2),
            "Req1": str(r1),
            "Req2": str(r2),
            "Req3": str(r3),
            "Req4": str(r4),
            "PT": str(z1),
            "ev_constraint": ev_constraint,
        }

        # self.debugging_data = {
        #     "constraint_balance": constrained,
        #     "constraint_problems": problems,
        #     "flow_balance": vflow,
        #     "flow_problems": vprob
        # }

    # ---------------------------
    # 5 — PRINT & STORE
    # ---------------------------

    def print_data(self, params, model_name):
        if params['ev_constraints']:
                
            def print_routes_ev():

                def print_vehicle_routes():

                    print("\n========================")
                    print("=== VEHICLE ROUTES ====")
                    print("========================\n")

                    # -----------------------------------------------------------
                    # 1. VEHICLE ROUTES
                    # -----------------------------------------------------------
                    for key, item in self.result.items():
                        if not key.startswith("Veh"):
                            continue

                        print(f"{key}:")

                        try:
                            route = ast.literal_eval(item)
                        except Exception as e:
                            print(f"  [Parse error for {key}]: {e}")
                            print()
                            continue

                        for step in route:
                            # Expecting:
                            # [ [i,j], i, node_type, req, T, soc, charge_time, pct ]
                            arc, node, node_type, req, T, soc, charge_time, pct = step

                            label = f"{node}"

                            if node_type == "pickup":
                                node_type_fmt = f"pickup({req})"
                            elif node_type == "dropoff":
                                node_type_fmt = f"dropoff({req})"
                            else:
                                node_type_fmt = node_type

                            base_line = f"  {label:<8} {node_type_fmt:<20} T={T:<7.1f} SOC={soc:5.1f}"

                            if node_type == "charging station":
                                base_line += f"   ChargeTime={charge_time:.1f}   +{pct:.1f}%"

                            print(base_line)

                        print()

                def print_request_routes():
                    # -----------------------------------------------------------
                    # 2. REQUEST ROUTES
                    # -----------------------------------------------------------
                    print("\n========================")
                    print("=== REQUEST ROUTES ====")
                    print("========================\n")

                    for key, item in self.result.items():
                        if not key.startswith("Req"):
                            continue

                        print(f"{key}:")

                        try:
                            req_list = ast.literal_eval(item)
                        except Exception as e:
                            print(f"  [Parse error for {key}]: {e}")
                            continue

                        # Format: [[r, [(i,m),(j,n)]], ...]
                        for segment in req_list:
                            if len(segment) != 2:
                                continue

                            r_id, arcs = segment
                            print(f" {arcs[0]} -> {arcs[1]}")
                        print()

                def print_PT_routes():
                    # -----------------------------------------------------------
                    # 3. PT ROUTES (Public Transport Arcs)
                    # -----------------------------------------------------------
                    print("\n========================")
                    print("=== PT ROUTES (z=1) ===")
                    print("========================\n")

                    if "PT" in self.result:
                        try:
                            pt_list = ast.literal_eval(self.result["PT"])
                        except Exception as e:
                            print(f"  [Parse error for PT]: {e}")
                            return

                        for elem in pt_list:
                            # elem = [((i,m),(j,n)), 'Departure X', 'T(i)=..', 'T(j)=..', 'z=1']
                            try:
                                arc, dep, Ti, Tj, z = elem
                                print(f"  {arc[0]} -> {arc[1]}   {dep},  {Ti},  {Tj},  {z}")
                            except:
                                print(f"  {elem}")

                    print("\n=== END OF ROUTES ===\n")

                def print_ev_constraints():

                    print("\n========================")
                    print("=== EV CONSTRAINTS ====")
                    print("========================\n")

                    if "ev_constraint" in self.result:

                        # If already a Python list, use it directly
                        if isinstance(self.result["ev_constraint"], list):
                            ev_list = self.result["ev_constraint"]
                        else:
                            # Fallback for string-stored lists
                            try:
                                ev_list = ast.literal_eval(self.result["ev_constraint"])
                            except Exception as e:
                                print(f"  [Parse error for ev_constraint]: {e}")
                                print("\n=== END OF EV CONSTRAINTS ===\n")
                                return


                        if not ev_list:
                            print("  (No active EV constraints)")
                            print("\n=== END OF EV CONSTRAINTS ===\n")
                            return

                        for elem in ev_list:

                            # Expected format:
                            # [
                            #   "ARC USED: (i,j)" or "ARC USED: vehicle k, (i,j)",
                            #   "  Battery at i: ...",
                            #   "  Battery at j: ...",
                            #   "  Theoretical consumption: ...",
                            #   "  Actual ΔB: ...",
                            #   "  Constraint: ..."
                            # ]

                            if not isinstance(elem, list):
                                print(f"  {elem} \n")
                                continue

                            for line in elem:
                                print(line)
                            print()

                    print("\n=== END OF EV CONSTRAINTS ===\n")
                
                print_vehicle_routes()
                print_request_routes()
                print_PT_routes()
                # print_ev_constraints()


            print_routes_ev()

        else:

            def print_routes_no_ev():

                print("\n========================")
                print("=== VEHICLE ROUTES ====")
                print("========================\n")

                # -----------------------------------------------------------
                # 1. VEHICLE ROUTES (no EV fields)
                # -----------------------------------------------------------
                for key, item in self.result.items():
                    if not key.startswith("Veh"):
                        continue

                    print(f"{key}:")

                    try:
                        route = ast.literal_eval(item)
                    except Exception as e:
                        print(f"  [Parse error for {key}]: {e}")
                        print()
                        continue

                    for step in route:
                        # EV OFF → format is:
                        # [ [i,j], i, node_type, req, T, soc ]
                        try:
                            arc, node, node_type, req, T, soc = step[:6]
                        except Exception:
                            print(f"  [Unexpected format in step: {step}]")
                            continue

                        label = f"{node}"

                        if node_type == "pickup":
                            node_type_fmt = f"pickup({req})"
                        elif node_type == "dropoff":
                            node_type_fmt = f"dropoff({req})"
                        else:
                            node_type_fmt = node_type

                        print(f"  {label:<8} {node_type_fmt:<20}  T={T:<7.1f} SOC={soc:5.1f}")

                    print()

                # -----------------------------------------------------------
                # 2. REQUEST ROUTES
                # -----------------------------------------------------------
                print("\n========================")
                print("=== REQUEST ROUTES ====")
                print("========================\n")

                for key, item in self.result.items():
                    if not key.startswith("Req"):
                        continue

                    print(f"{key}:")

                    try:
                        req_list = ast.literal_eval(item)
                    except Exception as e:
                        print(f"  [Parse error for {key}]: {e}")
                        continue

                    # [[req_id, [(i,m),(j,n)]], ...]
                    for segment in req_list:
                        if len(segment) != 2:
                            continue

                        r_id, arcs = segment
                        print(f"     {arcs[0]} -> {arcs[1]}")
                    print()

                # -----------------------------------------------------------
                # 3. PT ROUTES
                # -----------------------------------------------------------
                print("\n========================")
                print("=== PT ROUTES (z=1) ===")
                print("========================\n")

                if "PT" in self.result:
                    try:
                        pt_list = ast.literal_eval(self.result["PT"])
                    except Exception as e:
                        print(f"  [Parse error for PT]: {e}")
                        return

                    for elem in pt_list:
                        # elem = [((i,m),(j,n)), 'Departure x', Ti, Tj, 'z=1']
                        try:
                            arc, dep, Ti, Tj, z = elem
                            print(f"  {arc[0]} -> {arc[1]}   {dep},  {Ti},  {Tj},  {z}")
                        except:
                            print(f"  {elem}")

                print("\n=== END OF ROUTES ===\n")

            print_routes_no_ev()


        self.logger.log_dict(self.result, model_name)

        # streaming: write row immediately
        self._append_csv_row(self.result)

        extractor = DARPRouteExtractor(
            m=self.m, vars_=self.vars_, sets=self.sets, params=self.params, options=params
        )
        extractor.summarize()


        # print(self.debugging_data)

        self.results.append(self.result)
        # self.debugging_data_all.append(self.debugging_data)

    # ---------------------------
    # MAIN EXECUTION PIPELINES
    # ---------------------------

    def run_single_model(self, params, model_name="SingleModel"):
        self.build_model(params, model_name)

        # self.test_mip_start_feasibility(self.m)

        self.optimise_model(params, model_name)
        self.debug_model(params, model_name)
        self.print_data(params, model_name)

    def run_all_combinations(self):

        boolean_params = {
            "duplicate_transfers": [True, False],
            "arc_elimination": [True, False],
            "variable_substitution": [True, False],
            # "subtour_elimination": [True, False],
            "transfer_node_strengthening": [True, False],
            "ev_constraints": [True, False],
            "timetabled_departures": [True, False],
            "use_imjn": [True, False],
            "MoPS": [True, False]
        }

        combos = list(itertools.product(*boolean_params.values()))
        total = len(combos)
        pbar = tqdm(total=total, desc="Running Experiments", ncols=95)

        for idx, combo in enumerate(combos, start=1):
            params = dict(zip(boolean_params.keys(), combo))
            params.update({
            # "duplicate_transfers": [True, False],
            # "arc_elimination": [True, False],
            # "variable_substitution": [True, False],
            "subtour_elimination": False,
            # "transfer_node_strengthening": [True, False],
            # "ev_constraints": [True, False],
            # "timetabled_departures": [True, False],
            # "use_imjn": [True, False],
            # "MoPS": [True, False]
            })

            # (1) Timetabled PT only allowed with IMJN layer
            if params["timetabled_departures"] and not params["use_imjn"]:
                pbar.update(1)
                continue

            # (2) Duplicate transfers incompatible with IMJN
            if params["duplicate_transfers"] and params["use_imjn"]:
                pbar.update(1)
                continue

            # (3) Transfer strengthening incompatible with IMJN
            if params["transfer_node_strengthening"] and params["use_imjn"]:
                pbar.update(1)
                continue

            model_name = self._model_id(params)

            try:
                self.run_single_model(params, model_name)
            except Exception as e:
                self.logger.log(f"Error in {model_name}: {e}")
                error_row = {"Model": model_name, **params, "Objective": f"Error: {e}"}
                self._append_csv_row(error_row)
                self.results.append(error_row)

            pbar.update(1)

        pbar.close()
        self.logger.save_json()



# === Example usage ===
if __name__ == "__main__":
    TIME_LIMIT = 2 * 60 * 60
    bool_params_singular = {
        "duplicate_transfers": True,
        "arc_elimination": True,
        "variable_substitution": True,
        "subtour_elimination": False,
        "transfer_node_strengthening": False,
        "ev_constraints": False,
        "timetabled_departures": False,
        "use_imjn": False,
        "MoPS": False,
        
        "LNS": True
    }

    runner = DARPExperimentRunner(time_limit=TIME_LIMIT)
    runner.run_single_model(params=bool_params_singular)
    # runner.run_all_combinations()



