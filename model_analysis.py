import re

def analyze_model_report(input_file="test1_PT_timetabled.txt", output_file="model_analysis_timetabled.txt"):
    """
    Analyzes the Gurobi verification report and extracts:
    1. Variables x, v, y, z that are equal to 1
    2. All constraints involving those variables
    Writes the result into another text file.
    """
    with open(input_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    variables_equal_1 = []
    constraints_with_vars = []
    capture_constraints = False

    # Regex to capture variable names and values
    var_pattern = re.compile(r"^([xyvz]\[.*\])\s*=\s*([0-9eE\.\-\+]+)")

    # Step 1: Extract variables x, v, y, z = 1
    for line in lines:
        match = var_pattern.match(line.strip())
        if match:
            name, val = match.groups()
            try:
                if abs(float(val) - 1.0) < 1e-6:
                    variables_equal_1.append(name)
            except ValueError:
                continue

    # Step 2: Identify when constraints start
    for i, line in enumerate(lines):
        if "=== CONSTRAINTS ===" in line:
            capture_constraints = True
            constraints_section = lines[i+1:]
            break

    # Step 3: Find constraints that include any active variable
    current_constraint = []
    current_name = ""
    for line in constraints_section:
        if line.startswith("==="):
            continue
        elif line.startswith(">>> CONSTRAINT GROUP:"):
            continue
        elif re.match(r"^[A-Za-z0-9_∑\[\]]", line):
            # If we hit a new constraint name
            if current_constraint:
                full_text = "".join(current_constraint)
                # Check if it contains any variable = 1
                if any(var in full_text for var in variables_equal_1):
                    constraints_with_vars.append(full_text.strip())
                current_constraint = []
            current_constraint = [line]
            current_name = line.strip()
        else:
            current_constraint.append(line)

    # Append last one
    if current_constraint:
        full_text = "".join(current_constraint)
        if any(var in full_text for var in variables_equal_1):
            constraints_with_vars.append(full_text.strip())

    # Step 4: Write results to file
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("=== VARIABLES EQUAL TO 1 ===\n")
        for var in sorted(variables_equal_1):
            f.write(f"{var}\n")

        f.write("\n===================================================\n" * 3)
        f.write("=== CONSTRAINTS INVOLVING ACTIVE VARIABLES ===\n")
        for constr in constraints_with_vars:
            f.write("---------------------------------------------------\n")
            f.write(constr + "\n")

    print(f"✅ Analysis complete. Results written to '{output_file}'.")


# Run analysis
if __name__ == "__main__":
    analyze_model_report("C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes\\Stage-Inria-01-09-2025--28-02-2026\\test1_PT_timetabled.txt", "model_analysis_timetabled.txt")

# import re

def analyze_model_report_with_passenger_balance(
    input_file="C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes\\Stage-Inria-01-09-2025--28-02-2026\\test1_PT_timetabled.txt",
    output_file="model_analysis_passenger_balance_timetabled.txt"
):
    """
    Analyzes the Gurobi verification report and extracts:
    1. Variables x, v, y, z = 1
    2. Passenger balance constraints for visited nodes
    """
    with open(input_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # === STEP 1: Extract variables (x, v, y, z) ===
    var_pattern = re.compile(r"^([xyvz]\[.*\])\s*=\s*([0-9eE\.\-\+]+)")
    active_vars = []
    visited_nodes = set()

    for line in lines:
        match = var_pattern.match(line.strip())
        if match:
            name, val = match.groups()
            try:
                if abs(float(val) - 1.0) < 1e-6:
                    active_vars.append(name)
                    # Parse nodes for x and v variables
                    if name.startswith(("x[", "v[")):
                        node_pairs = re.findall(r"\(([^)]+)\)", name)
                        if node_pairs:
                            for node in node_pairs:
                                visited_nodes.add(node)
            except ValueError:
                continue

    # === STEP 2: Extract constraint section ===
    start_index = None
    for i, line in enumerate(lines):
        if "=== CONSTRAINTS ===" in line:
            start_index = i + 1
            break

    constraints = lines[start_index:] if start_index else []
    passenger_balance_constraints = []

    # === STEP 3: Find passenger balance constraints for visited nodes ===
    current_constr = []
    current_name = ""

    for line in constraints:
        # Skip separators and group headers
        if line.startswith("===") or line.startswith(">>>"):
            continue

        # Identify constraint name
        if re.match(r"^[∑A-Za-z0-9_]", line):
            # store previous constraint block
            if current_constr:
                full_text = "".join(current_constr)
                # Check if it's a passenger balance constraint and involves a visited node
                if "fi_r" in full_text and any(node in full_text for node in visited_nodes):
                    passenger_balance_constraints.append(full_text.strip())
                current_constr = []
            current_constr = [line]
            current_name = line.strip()
        else:
            current_constr.append(line)

    # handle the last constraint
    if current_constr:
        full_text = "".join(current_constr)
        if "fi_r" in full_text and any(node in full_text for node in visited_nodes):
            passenger_balance_constraints.append(full_text.strip())

    # === STEP 4: Write results ===
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("=== VARIABLES (x, v, y, z) EQUAL TO 1 ===\n")
        for v in sorted(active_vars):
            f.write(f"{v}\n")

        f.write("\n===================================================\n" * 3)
        f.write("=== VISITED NODES ===\n")
        for node in sorted(visited_nodes):
            f.write(f"{node}\n")

        f.write("\n===================================================\n" * 3)
        f.write("=== PASSENGER BALANCE CONSTRAINTS FOR VISITED NODES ===\n")

        for constr in passenger_balance_constraints:
            f.write("---------------------------------------------------\n")
            f.write(constr + "\n")

    print(f"✅ Analysis complete. Results written to '{output_file}'.")


# # Run analysis
# if __name__ == "__main__":
#     analyze_model_report_with_passenger_balance()
