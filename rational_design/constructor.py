import random
from pathlib import Path
from Bio import SeqIO

class LibraryConstructor:
    def construct(self, raw_target_dir, raw_bg_dir, output_paths, counts):
        print("--- [STAGE 1] BUILDING DATASETS (Strain-based) ---")
        for p in output_paths.values(): Path(p).parent.mkdir(parents=True, exist_ok=True)

        self._process_category("design", raw_target_dir, output_paths["design_target"], counts["design_target"])
        self._process_category("design", raw_bg_dir, output_paths["design_background"], counts["design_background"])
        self._process_category("validation", raw_target_dir, output_paths["validation_target"], counts["validation_target"])
        self._process_category("validation", raw_bg_dir, output_paths["validation_background"], counts["validation_background"])

    def _process_category(self, stage, input_dir, output_file, n_strains):
        path = Path(input_dir)
        if not path.exists():
            with open(output_file, "w") as f: pass
            return

        exts = ["*.fasta", "*.fa", "*.fna", "*.fas", "*.FNA", "*.FA", "*.FAS"]
        all_files = []
        if path.is_file():
            all_files = [path]
        else:
            for e in exts:
                all_files.extend(list(path.glob(e)))
        all_files = sorted(set(all_files))

        if not all_files:
            print(f"   ⚠️ No files found in {path.name}")
            with open(output_file, "w") as f: pass
            return

        print(f"   ➡ {stage.capitalize()} - {path.name}: Found {len(all_files)} files.")

        final_records = []
        for f in all_files:
            try:
                seqs = list(SeqIO.parse(f, "fasta"))
                for r in seqs:
                    r.id = f"{f.stem}|{r.id}"
                    final_records.append(r)
            except Exception as e:
                print(f"      ❌ Error parsing {f.name}: {e}")
                continue

            avg_len_mb = sum(len(r.seq) for r in seqs) / len(seqs) / 1_000_000 if seqs else 0
            print(f"      ✅ {f.stem} ({len(seqs)} contig(s), avg {avg_len_mb:.2f} MB): Loaded 1 strain.")

        limit = max(0, int(n_strains or 0))
        if limit > 0 and len(final_records) > limit:
            print(f"      📉 {stage.capitalize()} sampling: {limit} / {len(final_records)} records.")
            final_records = random.sample(final_records, limit)
        else:
            print(f"      📦 {stage.capitalize()} dataset: using all {len(final_records)} records.")

        with open(output_file, "w") as f: SeqIO.write(final_records, f, "fasta")
