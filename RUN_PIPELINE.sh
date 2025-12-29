#!/bin/bash
# Force UTF-8 Encoding for Emojis
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

clear
echo "========================================================"
echo "   🧬 RATIONAL PRIMER DESIGN PIPELINE (Linux / macOS) 🧬"
echo "========================================================"
echo

# --- STEP 1: PROJECT NAME ---
read -p "> Enter Project Name (e.g. MyTest): " proj_name
echo

# --- STEP 2: CONFIGURATION ---
echo " [CONFIGURATION]"
echo " -------------------------------------------------------"
echo " Default config: config/parameters.json"
echo " (Press Enter to use default, or drag a new JSON here)"
echo " -------------------------------------------------------"
read -p "> Drag Parameter JSON here: " param_file

# Clean drag-and-drop quotes (Handles Mac and Linux variations)
param_file="${param_file//\'/}"
param_file="${param_file//\"/}"

# Set default if empty
if [ -z "$param_file" ]; then
    param_file="config/parameters.json"
fi

echo
echo " [INFO] Using Configuration: $param_file"
echo

echo "========================================================"
echo " Please select your operation mode:"
echo
echo "   [1] FULL AUTOMATION (Download from NCBI)"
echo "   [2] LOCAL MODE (Use my own FASTA files)"
echo
echo "========================================================"
read -p "Enter choice (1 or 2): " mode

if [ "$mode" == "1" ]; then
    echo
    echo " [MODE: DOWNLOAD FROM NCBI]"
    echo " -------------------------------------------------------"
    echo " Tip: Drag and drop your JSON config files here."
    echo " -------------------------------------------------------"
    echo
    read -p "> Drag & Drop your TARGET JSON file here: " t_conf
    echo
    read -p "> Drag & Drop your BACKGROUND JSON file here: " b_conf
    
    # Clean quotes
    t_conf="${t_conf//\'/}"
    t_conf="${t_conf//\"/}"
    b_conf="${b_conf//\'/}"
    b_conf="${b_conf//\"/}"

    echo " 🚀 Starting Pipeline..."
    rational-design pipeline --target_config "$t_conf" --bg_config "$b_conf" --email "your_email@example.com" --out "$proj_name" --params "$param_file"

elif [ "$mode" == "2" ]; then
    echo
    echo " [MODE: LOCAL FILES]"
    echo " -------------------------------------------------------"
    echo " Default Target:     database/target"
    echo " Default Background: database/background"
    echo
    echo " (Press Enter to use defaults, or drag a specific folder)"
    echo " -------------------------------------------------------"
    echo

    # TARGET FOLDER
    read -p "> Path to TARGET folder: " t_path
    t_path="${t_path//\'/}"
    t_path="${t_path//\"/}"
    # Remove trailing slash if present
    t_path="${t_path%/}" 
    if [ -z "$t_path" ]; then t_path="database/target"; fi

    # BACKGROUND FOLDER
    read -p "> Path to BACKGROUND folder: " b_path
    b_path="${b_path//\'/}"
    b_path="${b_path//\"/}"
    b_path="${b_path%/}"
    if [ -z "$b_path" ]; then b_path="database/background"; fi

    echo
    echo " [INFO] Target Source:     $t_path"
    echo " [INFO] Background Source: $b_path"
    echo
    echo " 🚀 Starting Pipeline using local data..."
    rational-design pipeline --local_target "$t_path" --local_bg "$b_path" --out "$proj_name" --params "$param_file"

else
    echo "Invalid choice. Exiting."
    exit 1
fi

echo
echo "========================================================"
echo " ✅ DONE! Check the folder '$proj_name' for 'FINAL_ASSAY.csv'."
echo "========================================================"