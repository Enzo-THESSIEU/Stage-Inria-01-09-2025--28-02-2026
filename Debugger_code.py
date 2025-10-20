import pandas as pd
import json
import os

def save_model_debug_data(m, vars_, model_name="model_debug"):
    """
    Save all variable values and basic constraint information
    after optimizing a Gurobi model for troubleshooting.
    """

    # === 1️⃣ VARIABLES ===
    var_records = []
    for v in m.getVars():
        var_records.append({
            "Name": v.VarName,
            "Value": v.X if v.X is not None else None,
            "LB": v.LB,
            "UB": v.UB,
            "VType": v.VType,
            "RC": v.RC if hasattr(v, "RC") else None,  # reduced cost (if continuous)
        })
    df_vars = pd.DataFrame(var_records)
    df_vars.to_csv(f"{model_name}_vars.csv", index=False)
    print(f"✅ Saved {len(df_vars)} variables to {model_name}_vars.csv")

    # === 2️⃣ CONSTRAINTS ===
    constr_records = []
    for c in m.getConstrs():
        constr_records.append({
            "Name": c.ConstrName,
            "Sense": c.Sense,
            "RHS": c.RHS,
            "Slack": c.Slack if hasattr(c, "Slack") else None,
            "Pi": c.Pi if hasattr(c, "Pi") else None,  # dual value
        })
    df_constrs = pd.DataFrame(constr_records)
    df_constrs.to_csv(f"{model_name}_constrs.csv", index=False)
    print(f"✅ Saved {len(df_constrs)} constraints to {model_name}_constrs.csv")

    # === 3️⃣ PARAMETERS & SUMMARY ===
    summary = {
        "ModelName": m.ModelName,
        "ObjectiveValue": m.ObjVal if m.SolCount > 0 else None,
        "BestBound": m.ObjBound,
        "Status": m.Status,
        "NumVars": m.NumVars,
        "NumConstrs": m.NumConstrs,
        "NumNZs": m.NumNZs,
    }
    with open(f"{model_name}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"✅ Saved model summary to {model_name}_summary.json")

    # === 4️⃣ OPTIONALLY: Export .sol or .lp file ===
    m.write(f"{model_name}.sol")  # solution file (text)
    m.write(f"{model_name}.lp")   # LP format (for model inspection)

    print(f"✅ Full model snapshot saved for {model_name}")


def save_constraint_lhs_rhs(m, file_name="constraint_values.csv"):
    """
    Compute the evaluated LHS and RHS of each constraint in the Gurobi model.

    Creates a CSV file with:
        Constraint Name | LHS Value | RHS Value | Sense | Slack | Pi (dual)
    """

    data = []

    for constr in m.getConstrs():
        name = constr.ConstrName
        sense = constr.Sense
        rhs = constr.RHS

        # Evaluate LHS value numerically
        expr = m.getRow(constr)
        lhs_val = expr.getValue()  # substitute all variable values (from current solution)

        # Collect
        data.append({
            "Constraint": name,
            "LHS": lhs_val,
            "RHS": rhs,
            "Sense": sense,
            "Slack": getattr(constr, "Slack", None),
            "Pi": getattr(constr, "Pi", None)
        })

    df = pd.DataFrame(data)
    df.to_csv(file_name, index=False)
    print(f"✅ Saved constraint LHS/RHS evaluation to {file_name}")
    return df
