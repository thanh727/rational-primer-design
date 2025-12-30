import streamlit as st
import os
import sys
import json
import subprocess
import tkinter as tk
from tkinter import filedialog

# --- APP CONFIG ---
st.set_page_config(page_title="Rational Primer Design", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .reportview-container { margin-top: -2em; }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stDeployButton {display:none;}
</style>
""", unsafe_allow_html=True)

st.title("🧬 Rational Primer Design: Desktop App")
st.markdown("A high-performance pipeline for designing and validating TaqMan assays.")

# --- SIDEBAR: CONFIGURATION ---
with st.sidebar:
    st.header("⚙️ Advanced Configuration")
    
    # --- NEW: BLAST OPTION ADDED HERE ---
    st.subheader("Workflow Options")
    use_blast = st.checkbox("Enable Target Gene Annotation (BLAST)", value=True, 
                            help="If checked, the pipeline will blast candidates against the database to find gene names.")
    st.divider()
    # ------------------------------------

    with st.expander("🧬 Biological Parameters", expanded=True):
        min_sens = st.slider("Min Sensitivity (%)", 50.0, 100.0, 95.0, 0.1)
        min_cons = st.slider("Min Conservation (0-1)", 0.5, 1.0, 0.90, 0.01)
        max_xr = st.slider("Max Cross-Reactivity (%)", 0.0, 100.0, 5.0, 0.1)
        
        c1, c2 = st.columns(2)
        prod_min = c1.number_input("Min Product (bp)", value=100)
        prod_max = c2.number_input("Max Product (bp)", value=350)
        
        c3, c4 = st.columns(2)
        p_len = c3.number_input("Primer Length", value=20)
        max_mm = c4.number_input("Max Mismatches", value=2)

    with st.expander("💻 System & Sampling", expanded=False):
        cpu = st.number_input("CPU Cores (0=Auto)", value=0)
        st.markdown("**Sampling Sizes (0 = Use All)**")
        samp_dt = st.number_input("Design Target", value=0)
        samp_db = st.number_input("Design Background", value=100)
        samp_vt = st.number_input("Validate Target", value=0)
        samp_vb = st.number_input("Validate Background", value=200)

    st.markdown("---")
    st.info("Settings are saved automatically when you run the pipeline.")

# --- HELPER: FOLDER SELECTION CALLBACK ---
def select_folder_callback(session_key):
    """
    Opens a native OS dialog.
    IMPORTANT: Updates the specific session_state key used by the text_input.
    """
    try:
        # Create hidden root window
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1) # Make it float on top
        
        # Open dialog
        folder_path = filedialog.askdirectory(master=root)
        
        # Clean up
        root.destroy()
        
        # Update the key directly. Streamlit will re-render with this value.
        if folder_path:
            st.session_state[session_key] = folder_path
            
    except Exception as e:
        st.error(f"Error opening folder dialog: {e}")

# --- INITIALIZE SESSION STATE IF MISSING ---
if "path_t_val" not in st.session_state:
    st.session_state["path_t_val"] = ""
if "path_b_val" not in st.session_state:
    st.session_state["path_b_val"] = ""
if "l_out_val" not in st.session_state:
    st.session_state["l_out_val"] = os.path.join(os.getcwd(), "results_local")

# --- MAIN TABS ---
tab_auto, tab_local = st.tabs(["🤖 Automatic Mode (NCBI)", "📂 Local File Mode"])

# --- HELPER: RUN PIPELINE ---
def run_pipeline(cmd):
    """Runs the CLI command and streams output to Streamlit"""
    st.divider()
    st.header("📝 Execution Logs")
    
    log_area = st.empty()
    logs = []

    # --- ENCODING FIX FOR WINDOWS ---
    my_env = os.environ.copy()
    my_env["PYTHONIOENCODING"] = "utf-8"
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        text=True,
        bufsize=1,
        env=my_env,
        encoding="utf-8"
    )
    
    while True:
        try:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                clean_line = line.strip()
                logs.append(clean_line)
                log_area.code("\n".join(logs[-20:]), language="bash")
                print(clean_line) 
        except UnicodeDecodeError:
            continue
            
    if process.returncode == 0:
        st.success("✅ Pipeline Finished Successfully!")
        st.balloons()
    else:
        st.error("❌ Pipeline Failed. Check the logs above.")

# --- HELPER: SAVE CONFIG ---
def save_params():
    params = {
        "min_sensitivity": min_sens,
        "design_min_conservation": min_cons,
        "validation_max_cross_reactivity": max_xr,
        "product_size_min": prod_min,
        "product_size_max": prod_max,
        "primer_length": p_len,
        "max_mismatch": max_mm,
        "cpu_cores": cpu,
        "design_target_sampling_size": samp_dt,
        "design_background_sampling_size": samp_db,
        "validation_target_sampling_size": samp_vt,
        "validation_background_sampling_size": samp_vb,
        "design_max_candidates": 50,
        "enable_blast": use_blast  # <--- LINKED TO CHECKBOX
    }
    os.makedirs("config", exist_ok=True)
    path = "config/gui_params.json"
    with open(path, "w") as f:
        json.dump(params, f, indent=4)
    return path

# --- TAB 1: AUTO MODE ---
with tab_auto:
    st.subheader("Download Genomes from NCBI & Design")
    
    col1, col2 = st.columns(2)
    email = col1.text_input("NCBI Email (Required)", placeholder="email@example.com")
    project_name = col2.text_input("Project Name", value="Auto_Run_01")
    
    st.markdown("#### 🎯 Search Queries")
    t_query = st.text_input("Target Query", placeholder="e.g., Salmonella enterica[Org] AND complete genome")
    b_query = st.text_input("Background Query", placeholder="e.g., Escherichia coli[Org] AND complete genome")
    
    out_dir = st.text_input("Output Folder", value="results_auto")
    
    if st.button("🚀 Start Auto Pipeline", type="primary"):
        if not email or "@" not in email:
            st.error("Please enter a valid email.")
        elif not t_query or not b_query:
            st.error("Please enter both Target and Background queries.")
        else:
            param_file = save_params()
            
            t_conf = {"t1": [t_query, 0.0]}
            b_conf = {"b1": [b_query, 0.0]}
            
            os.makedirs("config", exist_ok=True)
            with open("config/t_conf.json", "w") as f: json.dump(t_conf, f)
            with open("config/b_conf.json", "w") as f: json.dump(b_conf, f)
            
            cmd = (f"\"{sys.executable}\" -u -m rational_design.cli pipeline "
                   f"--out \"{out_dir}\" "
                   f"--email \"{email}\" "
                   f"--target_config \"config/t_conf.json\" "
                   f"--bg_config \"config/b_conf.json\" "
                   f"--params \"{param_file}\"")
            
            run_pipeline(cmd)

# --- TAB 2: LOCAL MODE ---
with tab_local:
    st.subheader("Use Local FASTA Files")
    st.info("💡 Click 'Browse' to select folders. The path will appear automatically.")
    
    l_proj = st.text_input("Local Project Name", value="Local_Run_01")
    
    # --- TARGET FOLDER SELECTION ---
    col_t1, col_t2 = st.columns([4, 1])
    with col_t1:
        path_t = st.text_input("Path to Target Folder", key="path_t_val")
    with col_t2:
        st.write("") 
        st.write("") 
        st.button("📂 Browse", key="btn_t", on_click=select_folder_callback, args=("path_t_val",))

    # --- BACKGROUND FOLDER SELECTION ---
    col_b1, col_b2 = st.columns([4, 1])
    with col_b1:
        path_b = st.text_input("Path to Background Folder", key="path_b_val")
    with col_b2:
        st.write("") 
        st.write("") 
        st.button("📂 Browse", key="btn_b", on_click=select_folder_callback, args=("path_b_val",))
    
    # --- OUTPUT FOLDER SELECTION ---
    col_o1, col_o2 = st.columns([4, 1])
    with col_o1:
        l_out = st.text_input("Local Output Folder", key="l_out_val")
    with col_o2:
        st.write("") 
        st.write("") 
        st.button("📂 Browse", key="btn_o", on_click=select_folder_callback, args=("l_out_val",))
    
    st.divider()

    if st.button("🚀 Start Local Pipeline", type="primary"):
        if not path_t or not os.path.exists(path_t):
            st.error("❌ Target folder path is invalid or empty.")
        elif not path_b or not os.path.exists(path_b):
            st.error("❌ Background folder path is invalid or empty.")
        else:
            param_file = save_params()
            
            cmd = (f"\"{sys.executable}\" -u -m rational_design.cli pipeline "
                   f"--out \"{l_out}\" "
                   f"--local_target \"{path_t}\" "
                   f"--local_bg \"{path_b}\" "
                   f"--params \"{param_file}\"")
            
            run_pipeline(cmd)