# --- Launch Method ---
"To run this file don't run as usual but go inside the correct folder <below> then run <streamlit run Streamlit_app.py> in the terminal"
# cd "\Users\enzot\Documents\Césure\1ère césure inria Lille\Codes\Preliminary_model"
# streamlit run Streamlit_app.py

# --- Imports ---

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
# from Preliminary_model.Plotter_Pyvis import display_vehicle_routes_from_idarp
from Preliminary_model.Plotter_Plotly_express import display_vehicle_routes_from_idarp
import streamlit.components.v1 as components
import pandas as pd
import os

# --- Home Page ---
def home_page():
    st.title("The integrated dial-a-ride Problem")
    st.write("Welcome! Choose what you want to do:")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Run a singular Model"):
            st.session_state["page"] = "singular"
    with col2:
        if st.button("Run multiple Models"):
            st.session_state["page"] = "multiple"  # placeholder for later

import streamlit as st

# --- Singular Model Setup Page ---
def singular_model_page():
    st.header("Singular Model Setup")

    # --- Enhancement options ---
    categories = [
        "Duplicate Transfer Nodes",
        "Arc Elimination",
        "Variable Substitution",
        "Subtour Elimination",
        # "Transfer-node Strengthening",
        # "Electric Vehicles",
        "Timetabled Departure",
        "Use IMJN Nodes"
    ]

    # --- Objective functions ---
    functions = {
        r"f_1": r"\min \sum_{(i,m)\in\mathcal{N}} \sum_{(j,n)\in\mathcal{N}} c_{imjn} v_{imjn}",
        r"f_2": r"\min \sum_{r \in \mathcal{R}} \frac{(T_{(n+r,1)} - T_{(r,1)} - d_r)}{t_{r,n+r}}",
    }

    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("Select Enhancements")
        selections = {cat: st.checkbox(cat, value=False) for cat in categories}

        st.subheader("Define Objective Weights")

        variables = {}
        total_weight = 0.0

        # Input ω_i for each objective
        for i, (f_name, f_eq) in enumerate(functions.items(), start=1):
            st.markdown(f"**Objective {f_name}:**")
            st.latex(f_eq)
            ω = st.number_input(
                f"Weight ω{i} for {f_name}",
                min_value=0.0,
                max_value=1.0,
                value=0.0,
                step=0.05,
                key=f"var{i}_singular"
            )
            variables[f"ω{i}"] = ω
            total_weight += ω

        # Show combined weighted objective dynamically
        st.markdown("---")
        st.subheader("Weighted Sum Objective Function")

        latex_expr = " + ".join(
            [f"{variables[f'ω{i}']:.2f}\\,{f_name}" for i, f_name in enumerate(functions.keys(), start=1)]
        )
        st.latex(r"f = " + latex_expr)

        # Validate sum of weights
        if abs(total_weight - 1.0) > 1e-6:
            st.error(f"⚠️ The sum of weights must equal 1. Current sum = {total_weight:.2f}")
            can_run = False
        else:
            st.success("✅ The sum of weights equals 1. You can run the test.")
            can_run = True

    with col2:
        st.subheader("Action")

        run_button = st.button("Run Test", type="primary", disabled=not can_run)
        if run_button:
            st.session_state["singular_selections"] = selections
            st.session_state["singular_variables"] = variables
            st.session_state["page"] = "singular_results"
            st.success("✅ Configuration saved. Proceeding to results page...")

# --- Singular Model Results Page ---
def singular_results_page():
    st.header("Singular Model Results")

    selections = st.session_state.get("singular_selections", {})
    variables = st.session_state.get("singular_variables", {})

    # Build model name from selections (you can adjust naming to your pattern)
    model_name = "IDARP_Model_" + "_".join([
        f"{key[:3].lower()}{1 if val else 0}" for key, val in selections.items()
    ])

    csv_path = "C:\\Users\\enzot\\Documents\\Césure\\1ère césure inria Lille\\Codes\\results_full_2h_1.csv"

    # Check if model exists
    df = pd.read_csv(csv_path)
    if model_name in df["Model"].values:
        st.success(f"✅ Found existing simulation: **{model_name}**")
        html_file = f"{model_name.replace(' ', '_').lower()}_routes.html"

        # Generate if missing locally
        if not os.path.exists(html_file):
            display_vehicle_routes_from_idarp(csv_path, model_name)

        # Display the HTML network plot
        components.html(open(html_file, "r", encoding="utf-8").read(), height=750)
    else:
        st.warning(f"⚠️ No existing simulation found for {model_name}. You may run a new one.")

    if st.button("Back to Home"):
        st.session_state["page"] = "home"


# --- Multiple Models Setup Page ---
def multiple_model_page():
    st.header("Multiple Models Setup")

    num_models = st.slider("Number of Models", 1, 11, 3)

    categories = [
        "Duplicate Transfer Nodes",
        "Arc Elimination",
        "Variable Substitution",
        "Subtour Elimination",
        # "Transfer-node Strengthening",
        # "Electric Vehicles",
        "Timetabled Departure",
        "Use IMJN Nodes"
    ]

    if "model_offset" not in st.session_state:
        st.session_state["model_offset"] = 0

    max_visible = 3
    start_idx = st.session_state["model_offset"]
    end_idx = min(start_idx + max_visible, num_models)

    st.write(f"Configure categories and variables for **{num_models}** models (showing {start_idx+1} to {end_idx})")

    col_nav1, col_nav2, col_nav3 = st.columns([1, 8, 1])
    with col_nav1:
        if st.button("⬅", disabled=(start_idx == 0)):
            st.session_state["model_offset"] = max(0, start_idx - 1)
            st.rerun()
    with col_nav3:
        if st.button("➡", disabled=(end_idx >= num_models)):
            st.session_state["model_offset"] = min(num_models - max_visible, start_idx + 1)
            st.rerun()

    visible_models = range(start_idx, end_idx)
    cols = st.columns(len(visible_models))

    selections = st.session_state.get("multiple_selections", {})
    variables = st.session_state.get("multiple_variables", {})

    for i, col in zip(visible_models, cols):
        with col:
            st.subheader(f"Model {i+1}")
            if i not in selections:
                selections[i] = {}
            if i not in variables:
                variables[i] = {}

            for cat in categories:
                key = f"m{i}_{cat}"
                selections[i][cat] = st.checkbox(cat, key=key, value=selections[i].get(cat, False))

            st.write("Variables:")
            for v in range(1, 11):
                variables[i][f"var{v}"] = st.number_input(f"var{v}", value=0.0, key=f"m{i}_var{v}")

    st.session_state["multiple_selections"] = selections
    st.session_state["multiple_variables"] = variables

    st.markdown("---")
    if st.button("Run Tests", type="primary"):
        st.session_state["current_model"] = 0
        st.session_state["num_models"] = num_models
        st.session_state["page"] = "multiple_results"


# --- Multiple Models Results Page ---
def multiple_results_page():
    st.header("Multiple Models Results")

    selections = st.session_state.get("multiple_selections", {})
    variables = st.session_state.get("multiple_variables", {})
    current_model = st.session_state.get("current_model", 0)
    num_models = st.session_state.get("num_models", len(selections))

    st.subheader(f"Results for Model {current_model+1} of {num_models}")

    # --- Build model name (same pattern as singular) ---
    current_selection = selections.get(current_model, {})
    model_name = "IDARP_Model_" + "_".join([
        f"{key[:3].lower()}{1 if val else 0}" for key, val in current_selection.items()
    ])

    csv_path = r"C:\Users\enzot\Documents\Césure\1ère césure inria Lille\Codes\results_full.csv"

    # --- Check if model exists ---
    df = pd.read_csv(csv_path)
    if model_name in df["Model"].values:
        st.success(f"✅ Found existing simulation: **{model_name}**")
        html_file = f"{model_name.replace(' ', '_').lower()}_routes.html"

        # Generate if not present locally
        if not os.path.exists(html_file):
            display_vehicle_routes_from_idarp(csv_path, model_name)

        # Display the HTML visualization
        components.html(open(html_file, "r", encoding="utf-8").read(), height=750)
    else:
        st.warning(f"⚠️ No existing simulation found for {model_name}. You may run a new one.")

    # --- Navigation between models ---
    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        if st.button("⬅ Previous", disabled=(current_model == 0)):
            st.session_state["current_model"] -= 1
            st.rerun()
    with col3:
        if st.button("Next ➡", disabled=(current_model == num_models - 1)):
            st.session_state["current_model"] += 1
            st.rerun()

    st.markdown("---")
    if st.button("Back to Home"):
        st.session_state["page"] = "home"

        
# --- Main App ---
def main():
    if "page" not in st.session_state:
        st.session_state["page"] = "home"

    page = st.session_state["page"]

    if page == "home":
        home_page()
    elif page == "singular":
        singular_model_page()
    elif page == "singular_results":
        singular_results_page()
    elif page == "multiple":
        multiple_model_page()
    elif page == "multiple_results":
        multiple_results_page()


if __name__ == "__main__":
    main()
