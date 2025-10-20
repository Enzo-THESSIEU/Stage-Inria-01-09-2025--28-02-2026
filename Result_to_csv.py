import pandas as pd
import Runner as r

results = r.results

# Convert to list of dicts
rows = []
for name, (objval, bound, v1_route, v2_route, cpu_time, elapsed, excess_ride_time) in results.items():
    rows.append({
        "Model": name,
        "Objective": objval,
        "Lower Bound": bound,
        "Vehicle 1 Route": str(v1_route),
        "Vehicle 2 Route": str(v2_route),
        "CPU Time (s)": cpu_time,
        "Elapsed Time (s)": elapsed,
        "Total Excess User Ride Time": excess_ride_time
    })

df = pd.DataFrame(rows)
df.to_csv("results.csv", index=False)
print("Results saved to results.csv")
