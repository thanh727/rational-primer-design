import streamlit as st
import os
import sys
import json
import subprocess
import base64
import platform

# --- 1. APP CONFIG ---
st.set_page_config(page_title="Rational Primer Design", page_icon="🧬", layout="wide")

# --- TKINTER SAFETY CHECK ---
# Windows & Linux: Enable Browse Button (Try to import Tkinter)
# macOS (Darwin): Disable Browse Button (It crashes Streamlit)
SYSTEM_OS = platform.system()
HAS_TK = False

if SYSTEM_OS == "Darwin":
    # macOS: Force Disable to prevent crash
    HAS_TK = False
else:
    # Windows & Linux: Try to enable
    try:
        import tkinter as tk
        from tkinter import filedialog
        HAS_TK = True
    except ImportError:
        HAS_TK = False

# --- 2. INITIALIZE SESSION STATE ---
if "reset_id" not in st.session_state: st.session_state["reset_id"] = 0 

if "email_val" not in st.session_state: st.session_state["email_val"] = ""
if "path_t_val" not in st.session_state: st.session_state["path_t_val"] = ""
if "path_b_val" not in st.session_state: st.session_state["path_b_val"] = ""
if "l_out_val" not in st.session_state: st.session_state["l_out_val"] = os.path.join(os.getcwd(), "results_local")
if "out_dir_val" not in st.session_state: st.session_state["out_dir_val"] = "results_auto"
if "auto_proj_val" not in st.session_state: st.session_state["auto_proj_val"] = "Auto_Run_01"
if "local_proj_val" not in st.session_state: st.session_state["local_proj_val"] = "Local_Run_01"

if "target_list" not in st.session_state:
    st.session_state["target_list"] = [{"query": "", "size": 0.0, "count": 50}]
if "bg_list" not in st.session_state:
    st.session_state["bg_list"] = [{"query": "", "size": 0.0, "count": 100}]

# --- 3. HELPER FUNCTIONS ---

def reset_app():
    st.session_state["path_t_val"] = ""
    st.session_state["path_b_val"] = ""
    st.session_state["l_out_val"] = os.path.join(os.getcwd(), "results_local")
    st.session_state["out_dir_val"] = "results_auto"
    st.session_state["auto_proj_val"] = "Auto_Run_01"
    st.session_state["local_proj_val"] = "Local_Run_01"
    st.session_state["target_list"] = [{"query": "", "size": 0.0, "count": 50}]
    st.session_state["bg_list"] = [{"query": "", "size": 0.0, "count": 100}]
    st.session_state["reset_id"] += 1

def get_clickable_image_html(image_path, target_url):
    try:
        with open(image_path, "rb") as f:
            data = f.read()
        b64_data = base64.b64encode(data).decode()
        return f'<a href="{target_url}" target="_blank"><img src="data:image/png;base64,{b64_data}" style="width:100%; max-width:100%;"></a>'
    except Exception:
        return None

def select_folder_callback(session_key):
    if not HAS_TK:
        return
    try:
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder_path = filedialog.askdirectory(master=root)
        root.destroy()
        if folder_path:
            st.session_state[session_key] = folder_path
    except Exception as e:
        st.error(f"Error opening folder dialog: {e}")

# --- 4. CSS STYLING ---
st.markdown("""
<style>
    .reportview-container { margin-top: -2em; }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stDeployButton {display:none;}
    .stButton button { margin-top: 0px; }
    div.stCode > div > pre {
        max-height: 500px;
        overflow-y: auto !important;
        white-space: pre-wrap !important;
    }
    button[data-baseweb="tab"] div p {
        font-size: 1.2rem !important; 
        font-weight: 700 !important;
    }
    button[data-baseweb="tab"] {
        padding: 1rem 2rem !important;
        margin-right: 10px !important;
    }
</style>
""", unsafe_allow_html=True)

st.title("🧬 Rational Primer Design: Desktop App")
st.markdown("A high-performance pipeline for designing and validating TaqMan assays.")

# --- 5. SIDEBAR ---
with st.sidebar:
    if st.button("🗑️ Reset Search Fields", on_click=reset_app, help="Clear inputs (keeps Email) to start over."):
        st.rerun()
    
    st.divider()
    st.header("⚙️ Advanced Configuration")
    
    # --- TARGET MODE SELECTOR ---
    st.subheader("Target Biology")
    target_mode = st.radio(
        "Genomic Structure:",
        ("dsDNA (Bacteria/Euk)", "ssRNA (Viruses)"),
        help="dsDNA: Checks if Primer Sequence is in genome (Standard). ssRNA: Checks if Binding Site is in genome (Viral)."
    )
    
    st.subheader("Workflow Options")
    use_blast = st.checkbox("Enable Target Gene Annotation (BLAST)", value=True, 
                            help="If checked, the pipeline will blast candidates against the database to find gene names.")
    st.divider()

    with st.expander("🧬 Biological Parameters", expanded=True):
        min_sens = st.slider("Min Sensitivity (%)", 50.0, 100.0, 95.0, 0.1)
        min_cons = st.slider("Min Conservation (0-1)", 0.5, 1.0, 0.90, 0.01)
        
        # Tm Range
        primer_tm_range = st.slider(
            "Primer Melting Temp (°C)", 
            min_value=40.0, max_value=75.0, value=(50.0, 64.0), step=0.5,
            help="Select the minimum and maximum melting temperature for primers."
        )
        tm_min, tm_max = primer_tm_range

        max_xr = st.slider("Max Cross-Reactivity (%)", 0.0, 100.0, 5.0, 0.1)
        
        c1, c2 = st.columns(2)
        prod_min = c1.number_input("Min Product (bp)", value=100)
        prod_max = c2.number_input("Max Product (bp)", value=350)
        
        c3, c4 = st.columns(2)
        p_len = c3.number_input("Primer Length", value=20)
        max_mm = c4.number_input("Max Mismatches", value=2)

    with st.expander("💻 System & Sampling", expanded=False):
        # --- MODIFIED: MAX CANDIDATES SLIDER ---
        max_cand = st.slider(
            "Max Candidates to Test", 
            min_value=10, 
            max_value=2000, 
            value=500, 
            step=10,
            help="Higher values bridge larger gaps (good for viruses) but are slower. Lower values are faster (good for bacteria)."
        )
        # ---------------------------------------
        
        cpu = st.number_input("CPU Cores (0=Auto)", value=0)
        st.markdown("**Validation Sampling (0 = Use All)**")
        samp_dt = st.number_input("Design Target", value=0)
        samp_db = st.number_input("Design Background", value=100)
        samp_vt = st.number_input("Validate Target", value=0)
        samp_vb = st.number_input("Validate Background", value=200)

    st.markdown("---")
    st.info("Settings are saved automatically when you run the pipeline.")
    
    st.markdown("---")
    st.markdown("### 🌐 About")
    
    logo_file = None
    if os.path.exists("logo.png"): logo_file = "logo.png"
    elif os.path.exists("logo.jpg"): logo_file = "logo.jpg"
        
    if logo_file:
        html_code = get_clickable_image_html(logo_file, "https://genomessages.com/")
        if html_code:
            st.markdown(html_code, unsafe_allow_html=True)
        else:
            st.image(logo_file, use_container_width=True)

    st.markdown("**Developed by Thanh Nguyen, PhD**")
    st.markdown("Powered by **[https://genomessages.com](https://genomessages.com/)**")
    st.caption("Advanced Bioinformatics Solutions")

# --- 6. MAIN CONTENT ---
tab_auto, tab_local = st.tabs(["🤖 Automatic Mode (NCBI)", "📂 Local File Mode"])

def run_pipeline(cmd):
    st.divider()
    st.header("📝 Execution Logs")
    
    log_area = st.empty()
    logs = []

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
                full_log_text = "\n".join(logs[-5000:])
                log_area.code(full_log_text, language="bash")
        except UnicodeDecodeError:
            continue
            
    if process.returncode == 0:
        st.success("✅ Pipeline Finished Successfully!")
        st.balloons()
    else:
        st.error("❌ Pipeline Failed. Check the logs above.")

def save_params():
    # Convert UI selection to code param
    t_mode_code = "rna" if "ssRNA" in target_mode else "dna"
    
    params = {
        "target_mode": t_mode_code,
        "min_sensitivity": min_sens,
        "design_min_conservation": min_cons,
        "primer_tm_min": tm_min,
        "primer_tm_max": tm_max,
        "primer_opt_tm": (tm_min + tm_max) / 2,
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
        "design_max_candidates": max_cand, # <--- CONTROLLED BY SLIDER
        "enable_blast": use_blast
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
    email = col1.text_input("NCBI Email (Required)", placeholder="email@example.com", key="email_val")
    project_name = col2.text_input("Project Name", key="auto_proj_val")
    
    st.markdown("---")
    
    # --- DYNAMIC TARGET SECTION ---
    st.markdown("#### 🎯 Target Group (Inclusion)")
    for i, item in enumerate(st.session_state["target_list"]):
        c1, c2, c3, c4 = st.columns([5, 2, 2, 1])
        with c1:
            item["query"] = st.text_input(
                f"Target Query #{i+1}", 
                value=item["query"], 
                key=f"t_q_{i}_{st.session_state.reset_id}", 
                placeholder="e.g., Salmonella enterica[Org]..."
            )
        with c2:
            item["size"] = st.number_input(
                f"Min Size (Mb) #{i+1}", 
                value=item["size"], 
                min_value=0.0, 
                step=0.1, 
                key=f"t_s_{i}_{st.session_state.reset_id}"
            )
        # Count Column
        with c3:
            item["count"] = st.number_input(
                f"Max Count #{i+1}", 
                value=item.get("count", 50), 
                min_value=0, 
                step=1,
                key=f"t_c_{i}_{st.session_state.reset_id}",
                help="0 = Download All. 50 = Random 50."
            )
        with c4:
            st.write(""); st.write("")
            if i > 0: 
                if st.button("❌", key=f"del_t_{i}_{st.session_state.reset_id}"):
                    st.session_state["target_list"].pop(i)
                    st.rerun()

    if st.button("➕ Add Target Search", key="add_t"):
        st.session_state["target_list"].append({"query": "", "size": 0.0, "count": 50})
        st.rerun()

    st.markdown("---")

    # --- DYNAMIC BACKGROUND SECTION ---
    st.markdown("#### 🛡️ Background Group (Exclusion)")
    for i, item in enumerate(st.session_state["bg_list"]):
        c1, c2, c3, c4 = st.columns([5, 2, 2, 1])
        with c1:
            item["query"] = st.text_input(
                f"Background Query #{i+1}", 
                value=item["query"], 
                key=f"b_q_{i}_{st.session_state.reset_id}", 
                placeholder="e.g., Escherichia coli[Org]..."
            )
        with c2:
            item["size"] = st.number_input(
                f"Min Size (Mb) #{i+1}", 
                value=item["size"], 
                min_value=0.0, 
                step=0.1, 
                key=f"b_s_{i}_{st.session_state.reset_id}"
            )
        # Count Column
        with c3:
            item["count"] = st.number_input(
                f"Max Count #{i+1}", 
                value=item.get("count", 100), 
                min_value=0, 
                step=1,
                key=f"b_c_{i}_{st.session_state.reset_id}",
                help="0 = Download All."
            )
        with c4:
            st.write(""); st.write("")
            if i > 0:
                if st.button("❌", key=f"del_b_{i}_{st.session_state.reset_id}"):
                    st.session_state["bg_list"].pop(i)
                    st.rerun()

    if st.button("➕ Add Background Search", key="add_b"):
        st.session_state["bg_list"].append({"query": "", "size": 0.0, "count": 100})
        st.rerun()
    
    st.markdown("---")

    # Output & Run
    col_o1, col_o2 = st.columns([4, 1])
    with col_o1:
        out_dir = st.text_input("Output Folder", key="out_dir_val")
    with col_o2:
        st.write("") 
        st.write("") 
        if HAS_TK:
            st.button("📂 Browse", key="btn_auto_out", on_click=select_folder_callback, args=("out_dir_val",))
        else:
            st.empty() 
    
    if st.button("🚀 Start Auto Pipeline", type="primary"):
        t_valid = [x for x in st.session_state["target_list"] if x["query"].strip()]
        b_valid = [x for x in st.session_state["bg_list"] if x["query"].strip()]

        if not email or "@" not in email:
            st.error("Please enter a valid email.")
        elif not t_valid:
            st.error("❌ You must provide at least one valid Target query.")
        elif not b_valid:
            st.error("❌ You must provide at least one valid Background query.")
        else:
            param_file = save_params()
            
            t_conf = {}
            for idx, item in enumerate(t_valid):
                # Pass 3 items now: [Query, Size, Count]
                t_conf[f"t{idx+1}"] = [item["query"], item["size"], item["count"]]
            
            b_conf = {}
            for idx, item in enumerate(b_valid):
                b_conf[f"b{idx+1}"] = [item["query"], item["size"], item["count"]]
            
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
    if HAS_TK:
        st.info("💡 Click 'Browse' to select folders.")
    else:
        st.info("💡 **Mac Users:** You can drag & drop the folder into the box below (if supported) or copy the path.")
    
    l_proj = st.text_input("Local Project Name", key="local_proj_val")
    
    col_t1, col_t2 = st.columns([4, 1])
    with col_t1:
        path_t = st.text_input("Path to Target Folder", key="path_t_val")
    with col_t2:
        st.write("") 
        st.write("") 
        if HAS_TK:
            st.button("📂 Browse", key="btn_t", on_click=select_folder_callback, args=("path_t_val",))

    col_b1, col_b2 = st.columns([4, 1])
    with col_b1:
        path_b = st.text_input("Path to Background Folder", key="path_b_val")
    with col_b2:
        st.write("") 
        st.write("") 
        if HAS_TK:
            st.button("📂 Browse", key="btn_b", on_click=select_folder_callback, args=("path_b_val",))
    
    col_o1, col_o2 = st.columns([4, 1])
    with col_o1:
        l_out = st.text_input("Local Output Folder", key="l_out_val")
    with col_o2:
        st.write("") 
        st.write("") 
        if HAS_TK:
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