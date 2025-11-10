##### Imports #####
import gurobipy as gb


class DARPRouteExtractor_HallPosada:
    """
    Extracts and reconstructs vehicle, request, and PT routes
    from the Häll & Posada (2017) Model 2 formulation.
    Compatible with:
      - homogeneous fleet
      - single resource
      - timetabled PT system
    """

    def __init__(self, m, vars_, sets, params):
        self.m = m
        self.vars_ = vars_
        self.sets = sets
        self.params = params

    # ----------------------------------------------------------------------
    def base(self, x):
        """Return flat node id from (i,m) or int."""
        return x[0] if isinstance(x, tuple) else x

    # ----------------------------------------------------------------------
    def extract_vehicle_routes(self):
        """
        Extract vehicle routes from x[i,m,j,n,k].
        Returns dict {k: [(i,m)->(j,n), ...]}.
        """
        x = self.vars_.get("x", {})
        K = self.sets["K"]
        A = self.sets["A"]
        zeroDepot = self.sets["zeroDepot"]
        endDepot = self.sets["endDepot"]

        routes = {}
        for k in K:
            arcs_used = [(i, j) for (i, j, kk) in x.keys()
                         if kk == k and x[i, j, kk].X > 0.5]
            if not arcs_used:
                routes[k] = []
                continue

            # Start at depot
            start_arc = next(((i, j) for (i, j) in arcs_used
                              if self.base(i) == self.base(zeroDepot[0])), None)
            if not start_arc:
                print(f"[Warning] Vehicle {k} has no outgoing arc from depot.")
                routes[k] = []
                continue

            route = [start_arc]
            current = start_arc[1]
            arcs_used.remove(start_arc)

            stuck = 0
            while self.base(current) != self.base(endDepot[0]) and stuck < 200:
                found_next = False
                for (i, j) in list(arcs_used):
                    if i == current:
                        route.append((i, j))
                        current = j
                        arcs_used.remove((i, j))
                        found_next = True
                        break
                if not found_next:
                    stuck += 1
                    if stuck > 10:
                        print(f"⚠️ Vehicle {k} route interrupted at node {current}")
                        break

            routes[k] = route
        return routes

    # ----------------------------------------------------------------------
    def extract_request_routes(self):
        """
        Extract individual passenger routes from y[r,i,m,j,n].
        Returns dict {r: [(i,m)->(j,n), ...]}.
        """
        y = self.vars_.get("y", {})
        t = self.vars_.get("t", {})
        R = self.sets["R"]

        routes = {r: [] for r in R}
        for (r, i, j) in y.keys():
            if y[r, i, j].X > 0.5:
                routes[r].append((i, j))

        # Sort arcs by start service time
        for r in R:
            routes[r].sort(key=lambda arc: t[arc[0]].X if arc[0] in t else 0)
        return routes

    # ----------------------------------------------------------------------
    def extract_PT_routes(self):
        """
        Extract all used Public Transport arcs z[i,j,d,r].
        Returns list of tuples (r, (i,m)->(j,n), departure, arrival).
        """
        z = self.vars_.get("z", {})
        TD, TA = self.params.get("TD", {}), self.params.get("TA", {})

        used_PT = []
        for (i, j, d, r) in z.keys():
            if z[i, j, d, r].X > 0.5:
                dep = TD.get((self.base(i), self.base(j), d))
                arr = TA.get((self.base(i), self.base(j), d))
                used_PT.append((r, (i, j), d, dep, arr))
        return used_PT

    # ----------------------------------------------------------------------
    def summarize(self):
        """Prints readable summaries of extracted routes."""
        veh_routes = self.extract_vehicle_routes()
        req_routes = self.extract_request_routes()
        pt_routes = self.extract_PT_routes()

        print("\n=== VEHICLE ROUTES ===")
        for k, route in veh_routes.items():
            if route:
                chain = " → ".join([f"{i}->{j}" for (i, j) in route])
                print(f"Vehicle {k}: {chain}")
            else:
                print(f"Vehicle {k}: [no route]")

        print("\n=== REQUEST ROUTES ===")
        for r, arcs in req_routes.items():
            if arcs:
                print(f"Request {r}: " + " → ".join([f"{i}->{j}" for (i, j) in arcs]))
            else:
                print(f"Request {r}: [no assigned path]")

        print("\n=== PUBLIC TRANSPORT ARCS USED ===")
        if not pt_routes:
            print("No PT arcs used.")
        else:
            for (r, (i, j), d, dep, arr) in pt_routes:
                print(f"Request {r}: {i}→{j} (dep={dep}, arr={arr})")

        return veh_routes, req_routes, pt_routes


# ----------------------------------------------------------------------
# Example usage (after solving the model)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    from Model_Runner_Hall_Posada import DARPModelBuilder_HallPosada

    builder = DARPModelBuilder_HallPosada(model_name="Hall_Posada_ExtractTest", TIME_LIMIT=5 * 60)
    model, vars_ = builder.build()
    builder.run()

    extractor = DARPRouteExtractor_HallPosada(model, vars_, builder.sets, builder.params)
    veh_routes, req_routes, pt_routes = extractor.summarize()
