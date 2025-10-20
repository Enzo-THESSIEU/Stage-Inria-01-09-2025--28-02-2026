import pandas as pd
import ast

# === USER INPUT ===
CSV_FILE = "results_full.csv"

# --- (optional) known shortest-path travel times per request (you can adjust manually) ---
# If you don't know the theoretical min travel times, set to 0 (it will just give ride duration)
shortest_times = {
    1: 100,   # request 1
    2: 120,   # request 2
    3: 130,   # request 3
    4: 150,   # request 4
}

def parse_route(route_str):
    """Safely convert a route string like '[[[0,1],0,...]]' to a Python list."""
    try:
        return ast.literal_eval(route_str)
    except Exception:
        return []

def get_pickup_and_dropoff_times(route):
    """Extract pickup and dropoff times per request from a parsed route."""
    pickups = {}
    dropoffs = {}

    for step in route:
        if len(step) < 5:
            continue
        _, _, action, req_id, time = step
        if action == "pickup":
            pickups[req_id] = time
        elif action == "dropoff":
            dropoffs[req_id] = time

    return pickups, dropoffs

def compute_excess_times(df):
    """Compute excess ride times for all rows in the dataframe."""
    results = []

    for idx, row in df.iterrows():
        model_name = row["Model"]

        # parse both vehicle routes
        route1 = parse_route(row["Vehicle 1 Route"])
        route2 = parse_route(row["Vehicle 2 Route"])

        # merge pickups/dropoffs from both
        p1, d1 = get_pickup_and_dropoff_times(route1)
        p2, d2 = get_pickup_and_dropoff_times(route2)
        pickups = {**p1, **p2}
        dropoffs = {**d1, **d2}

        # compute ride and excess times
        for req_id in sorted(set(pickups) & set(dropoffs)):
            ride_time = dropoffs[req_id] - pickups[req_id]
            excess = ride_time - shortest_times.get(req_id, 0)
            results.append({
                "Model": model_name,
                "Request": req_id,
                "Pickup": pickups[req_id],
                "Dropoff": dropoffs[req_id],
                "Ride Time": round(ride_time, 2),
                "Excess Ride Time": round(excess, 2),
            })

    return pd.DataFrame(results)

# === RUN ===
df = pd.read_csv(CSV_FILE)
excess_df = compute_excess_times(df)

# === OUTPUT ===
print("\n=== Excess User Ride Times ===")
print(excess_df)

# Save to CSV
excess_df.to_csv("excess_user_ride_times.csv", index=False)
print("\n✅ Saved to excess_user_ride_times.csv")
