import random
import re
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

    def _collect_fasta_files(self, input_dir):
        path = Path(input_dir)
        exts = ["*.fasta", "*.fa", "*.fna", "*.fas", "*.FNA", "*.FA", "*.FAS"]
        if path.is_file():
            return [path]
        files = []
        for e in exts:
            files.extend(list(path.glob(e)))
        return sorted(set(files))

    def _parse_records(self, file_path):
        seqs = list(SeqIO.parse(file_path, "fasta"))
        for r in seqs:
            r.id = f"{file_path.stem}|{r.id}"
        return seqs

    def _group_key(self, filename):
        prefix = filename.stem.split("_")[0]
        if re.match(r"^[bt]\d+$", prefix):
            return prefix
        return None

    def _process_category(self, stage, input_dir, output_file, n_strains):
        path = Path(input_dir)
        if not path.exists():
            with open(output_file, "w") as f: pass
            return

        all_files = self._collect_fasta_files(path)
        if not all_files:
            print(f"   ⚠️ No files found in {path.name}")
            with open(output_file, "w") as f: pass
            return

        # Determine grouping: NCBI mode (b1,b2,...) → per-group; local mode → single pool
        group_keys = {self._group_key(f) for f in all_files}
        has_groups = any(k is not None for k in group_keys)

        if has_groups:
            groups: dict[str, list[Path]] = {}
            for f in all_files:
                key = self._group_key(f)
                groups.setdefault(key, []).append(f)
        else:
            groups = {"_all": all_files}

        limit = max(0, int(n_strains or 0))
        final_records = []

        for grp_name, grp_files in sorted(groups.items()):
            file_records: list[tuple[str, list]] = []
            for f in grp_files:
                try:
                    seqs = self._parse_records(f)
                    file_records.append((f.stem, seqs))
                except Exception as e:
                    print(f"   ❌ Error parsing {f.name}: {e}")
                    continue

            n_files = len(file_records)
            sampled = limit > 0 and n_files > limit

            if stage == "design" and "non" in path.name.lower():
                strategy = "sampled background for discovery" if sampled else "full background for discovery"
            elif stage == "design":
                strategy = "full target panel for discovery"
            elif stage == "validation" and "non" in path.name.lower():
                strategy = "sampled background for final specificity" if sampled else "full background panel for final specificity"
            else:
                strategy = "full target panel for final validation"

            n_keep = min(limit, n_files) if sampled else n_files
            print(f"\n   {path.name} / {grp_name}:")
            print(f"      strategy: {strategy}")
            print(f"      selected: {n_keep} / {n_files} strains")

            if sampled:
                selected = random.sample(file_records, limit)
                for stem, seqs in selected:
                    final_records.extend(seqs)
            else:
                for stem, seqs in file_records:
                    final_records.extend(seqs)

        print(f"\n   📊 {path.name}: {len(final_records)} contig(s) from final strains.\n")

        with open(output_file, "w") as f:
            SeqIO.write(final_records, f, "fasta")
