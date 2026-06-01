import json
import os
import subprocess
from datetime import datetime

payload = {
  "action": "propose_multiplex",
  "run_immediately": True,
  "targets": [
    {
      "query": "Streptococcus pneumoniae",
      "count": 500,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Haemophilus influenzae",
      "count": 500,
      "size": 0.0,
      "type": "genome"
    }
  ],
  "background": [
    {
      "query": "Streptococcus pyogenes",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Streptococcus suis",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Haemophilus parainfluenzae",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Neisseria meningitidis",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Moraxella catarrhalis",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Klebsiella pneumoniae",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    },
    {
      "query": "Staphylococcus aureus",
      "count": 50,
      "size": 0.0,
      "type": "genome"
    }
  ]
}

def run_multiplex(payload):
    email = "thanh.nt@v-extreme.com" # Dummy email for NCBI
    bg_dict = {}
    for i, bg in enumerate(payload["background"]):
        bg_dict[f"bg_{i}"] = [bg["query"], bg["size"], bg["count"], bg.get("type", "genome")]
    
    bg_json = "/Users/thanhnt/Desktop/Primer_tools/rational_primer_design/runs/bg_config.json"
    with open(bg_json, "w") as f:
        json.dump(bg_dict, f, indent=4)
        
    designed_folders = []
    
    for t_idx, target in enumerate(payload["targets"]):
        t_query = target["query"]
        t_dict = {
            "target": [t_query, target["size"], target["count"], target.get("type", "genome")]
        }
        
        t_json = "/Users/thanhnt/Desktop/Primer_tools/rational_primer_design/runs/target_config.json"
        with open(t_json, "w") as f:
            json.dump(t_dict, f, indent=4)
            
        # Create output directory
        safe_name = t_query.replace(" ", "_").replace('"', '').replace('[Organism]', '').replace('AND', '').strip()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = f"/Users/thanhnt/Desktop/Primer_tools/rational_primer_design/runs/{safe_name}/run_{timestamp}"
        
        cmd = [
            "python", "-m", "rational_design.cli", "pipeline",
            "--out", out_dir,
            "--email", email,
            "--target_config", t_json,
            "--bg_config", bg_json,
            "--ai_base_url", "http://localhost:11434/v1",
            "--ai_model", "llama3",
            "--language", "en"
        ]
        
        print(f"Running pipeline for target {t_idx+1}/{len(payload['targets'])}: {t_query}")
        print(" ".join(cmd))
        
        subprocess.run(cmd)
        designed_folders.append(out_dir)
        
    # Run multiplex analyze
    if len(designed_folders) > 1:
        multi_out = "/Users/thanhnt/Desktop/Primer_tools/rational_primer_design/runs/Multiplex_Result_" + datetime.now().strftime("%Y%m%d_%H%M%S")
        os.makedirs(multi_out, exist_ok=True)
        multi_cmd = [
            "python", "-m", "rational_design.cli", "multiplex_analyze",
            "-f"
        ] + designed_folders + [
            "-o", multi_out,
            "--assay_type", "qPCR",
            "--ai_base_url", "http://localhost:11434/v1",
            "--ai_model", "llama3",
            "--language", "en"
        ]
        print(f"Running Multiplex Analysis: {' '.join(multi_cmd)}")
        subprocess.run(multi_cmd)
        print(f"Multiplex Result saved in {multi_out}")

if __name__ == "__main__":
    os.makedirs("/Users/thanhnt/Desktop/Primer_tools/rational_primer_design/runs", exist_ok=True)
    run_multiplex(payload)
