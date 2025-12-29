import os
import shutil
import random
from pathlib import Path
from Bio import SeqIO

class LibraryConstructor:
    def __init__(self):
        pass

    def construct(self, raw_target_dir, raw_bg_dir, output_paths, counts):
        print("==========================================")
        print("🧬 PHASE 1: GENERATING DATASETS")
        print("   (Validation data now includes Design data)")
        print("==========================================")

        for p in output_paths.values():
            Path(p).parent.mkdir(parents=True, exist_ok=True)

        self._process_category("design", raw_target_dir, output_paths["design_target"], counts["design_target"], {})
        self._process_category("design", raw_bg_dir, output_paths["design_background"], counts["design_background"], {})

        print("\n==========================================")
        print("🧪 PHASE 2: GENERATING VALIDATION DATASETS")
        print("==========================================")
        
        self._process_category("validation", raw_target_dir, output_paths["validation_target"], counts["validation_target"], {})
        self._process_category("validation", raw_bg_dir, output_paths["validation_background"], counts["validation_background"], {})

        print("\n✨ Library Construction Complete.")

    def _process_category(self, stage_name, input_dir, output_file, n_sample, exclude_indices):
        print(f"\n📂 Processing Folder: {Path(input_dir).name}")
        print(f"   ➡ {stage_name.capitalize()} Output: {Path(output_file).name}")
        print(f"   ➡ Request: {n_sample if n_sample > 0 else 'ALL'} seqs/file")

        input_path = Path(input_dir)
        if not input_path.exists():
            print(f"❌ Error: Directory not found: {input_path}")
            return {}

        all_records = []
        new_usage = {}
        files = list(input_path.glob("*.fasta")) + list(input_path.glob("*.fa"))
        
        for fasta_file in files:
            try:
                recs = list(SeqIO.parse(fasta_file, "fasta"))
                total_in_file = len(recs)
                available_indices = list(range(total_in_file))
                
                if n_sample <= 0 or n_sample >= len(available_indices):
                    chosen_indices = available_indices
                    msg = f"All {len(chosen_indices)}"
                else:
                    chosen_indices = random.sample(available_indices, n_sample)
                    msg = f"Sampled {n_sample}/{len(available_indices)}"

                for idx in chosen_indices:
                    r = recs[idx]
                    r.id = f"{r.id}|{fasta_file.stem}" 
                    r.description = ""
                    all_records.append(r)
                
                new_usage[fasta_file.name] = set(chosen_indices)
                print(f"   ✔ {fasta_file.name}: {msg}")

            except Exception as e:
                print(f"   ⚠️ Skipping {fasta_file.name}: {e}")

        if all_records:
            with open(output_file, "w") as f:
                SeqIO.write(all_records, f, "fasta")
            print(f"   ✅ Written {len(all_records)} total sequences to {Path(output_file).name}")
        else:
            print(f"   ⚠️ No sequences found for {stage_name}.")
            with open(output_file, "w") as f: pass

        return new_usage