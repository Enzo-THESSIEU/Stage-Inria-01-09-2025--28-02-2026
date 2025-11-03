def write_model_verification_report(m, file_path="model_verification.txt"):
    """
    Write a complete model verification report listing every variable
    and constraint with their values after optimization.
    Displays five-line separators between variable/constraint groups.
    """
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"Model: {m.ModelName}\n")
        f.write(f"Objective Value: {m.ObjVal:.6f}\n\n")

        # === VARIABLES ===
        f.write("=== VARIABLES ===\n")
        prev_prefix = None

        for var in sorted(m.getVars(), key=lambda v: v.VarName):
            prefix = var.VarName.split("[", 1)[0] if "[" in var.VarName else var.VarName

            # Separator when changing variable group
            if prefix != prev_prefix:
                if prev_prefix is not None:
                    f.write("===================================================\n" * 5)
                f.write(f"\n>>> VARIABLE GROUP: {prefix}\n")
                prev_prefix = prefix

            try:
                f.write(f"{var.VarName} = {var.X:.6f}\n")
            except Exception:
                f.write(f"{var.VarName} = (unavailable)\n")

        f.write("\n===================================================\n" * 5)
        f.write("=== CONSTRAINTS ===\n")

        # === CONSTRAINTS ===
        prev_prefix = None
        for constr in sorted(m.getConstrs(), key=lambda c: c.ConstrName):
            prefix = constr.ConstrName.split("[", 1)[0] if "[" in constr.ConstrName else constr.ConstrName

            # Separator when constraint group changes
            if prefix != prev_prefix:
                if prev_prefix is not None:
                    f.write("===================================================\n" * 5)
                f.write(f"\n>>> CONSTRAINT GROUP: {prefix}\n")
                prev_prefix = prefix

            try:
                slack = constr.Slack
                f.write(
                    "===================================================\n"
                    f"{constr.ConstrName}:\n"
                    f"LHS = {m.getRow(constr)} "
                    f"→ Value = {m.getRow(constr).getValue():.6f}, "
                    f"Sense = {constr.Sense}, "
                    f"RHS = {constr.RHS:.6f}, "
                    f"Slack = {slack:.6f}\n"
                )
            except Exception:
                f.write(f"{constr.ConstrName}: (error reading constraint)\n")

        f.write("\n===================================================\n" * 5)
        f.write("=== END OF REPORT ===\n")

    print(f"✅ Model verification report written to {file_path}")
