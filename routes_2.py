def base(node):
    """Return the base node ID (for IMJN tuples)."""
    return node[0] if isinstance(node, tuple) else node


def safe_node_info(nodes, i):
    """Safely fetch node info (type, request id) for node i."""
    i_base = base(i)
    try:
        row = nodes[nodes[:, 0] == i_base][0]
        node_type = row[1]
        req_id = row[2] if len(row) > 2 else ""
    except Exception:
        node_type, req_id = "unknown", ""
    return node_type, req_id


def extract_route_final(vars_, nodes, N, K, variable_substitution=True, use_imjn=False):
    """
    Universal route extractor — works with:
    - variable_substitution: True (uses v) or False (uses x)
    - use_imjn: True (artificial node tuples) or False (flat node IDs)
    Orders request routes by service time (T_node).
    """

    y = vars_["y"]
    T_node = vars_["T_node"]

    # === VEHICLE ARCS ===
    used_arcs = {}

    if variable_substitution:
        v = vars_["v"]
        for (i, j) in v.keys():
            if v[(i, j)].X > 0.5:
                node_type, req_id = safe_node_info(nodes, i)
                used_arcs[(i, j)] = [
                    [i, j],
                    i,
                    node_type,
                    req_id,
                    T_node[j].X if (j in T_node) else None,
                ]
    else:
        x = vars_["x"]
        for (k, i, j) in x.keys():
            if x[(k, i, j)].X > 0.5:
                node_type, req_id = safe_node_info(nodes, i)
                used_arcs[(i, j)] = [
                    [i, j],
                    i,
                    node_type,
                    req_id,
                    T_node[j].X if (j in T_node) else None,
                ]

    # --- Detect first trips (from depot) ---
    used_arcs_vehicle_1, used_arcs_vehicle_2 = [], []
    first_trips = [used_arcs.pop(k) for k in list(used_arcs.keys()) if base(k[0]) == 0]

    if first_trips:
        used_arcs_vehicle_1.append(first_trips[0])
    if len(first_trips) > 1:
        used_arcs_vehicle_2.append(first_trips[1])

    # --- Stitch vehicle paths ---
    stuck_counter = 0
    while used_arcs and stuck_counter < 100:
        progress = False
        for key in list(used_arcs.keys()):
            i, j = key
            if used_arcs_vehicle_1 and i == used_arcs_vehicle_1[-1][0][1]:
                used_arcs_vehicle_1.append(used_arcs.pop(key))
                progress = True
            elif used_arcs_vehicle_2 and i == used_arcs_vehicle_2[-1][0][1]:
                used_arcs_vehicle_2.append(used_arcs.pop(key))
                progress = True
        if not progress:
            stuck_counter += 1
            if stuck_counter == 3:
                print("⚠️ No progress in vehicle stitching — breaking loop.")
                break

    # === REQUEST ARCS ===
    used_arcs_y = {
        (int(r), i, j): [int(r), [i, j]]
        for (r, i, j) in y.keys()
        if y[(r, i, j)].X > 0.5
    }

    # Group arcs per request
    req_routes = {r: [] for r in set(r for (r, _, _) in used_arcs_y.keys())}

    # Identify earliest node for each request
    for (r, i, j), val in list(used_arcs_y.items()):
        if not req_routes[r] or (
            T_node[i].X if i in T_node else float("inf")
        ) < (T_node[req_routes[r][0][1][0]].X if req_routes[r] else float("inf")):
            req_routes[r].insert(0, val)
            used_arcs_y.pop((r, i, j), None)

    # --- Stitch sequentially ---
    iter_limit = 200
    for _ in range(iter_limit):
        if not used_arcs_y:
            break
        progress = False
        for key in list(used_arcs_y.keys()):
            r, i, j = key
            if req_routes[r]:
                last_j = req_routes[r][-1][1][1]
                if i == last_j:
                    req_routes[r].append(used_arcs_y.pop(key))
                    progress = True
        if not progress:
            break

    # Append leftover arcs (if disjoint)
    if used_arcs_y:
        print(f"⚠️ Unconnected arcs left for some requests: {len(used_arcs_y)}")
        for key in used_arcs_y.keys():
            r, i, j = key
            req_routes.setdefault(r, []).append(used_arcs_y[key])

    # === SORT REQUEST ROUTES BY SERVICE TIME ===
    for r in req_routes:
        req_routes[r].sort(
            key=lambda a: T_node[a[1][1]].X if a[1][1] in T_node else float("inf")
        )

    # === Return ordered routes ===
    return (
        used_arcs_vehicle_1,
        used_arcs_vehicle_2,
        req_routes.get(1, []),
        req_routes.get(2, []),
        req_routes.get(3, []),
        req_routes.get(4, []),
    )
