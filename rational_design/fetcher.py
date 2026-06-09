import os
import time
import random
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

class SequenceFetcher:
    """
    A robust NCBI downloader that handles large datasets, network retries,
    dynamic genome size filtering, AND RANDOM SAMPLING.
    """
    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        chunk_size: int = 200,
        max_retries: int = 3
    ):
        self.email = email
        self.chunk_size = chunk_size
        self.max_retries = max_retries

        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def fetch_accession_numbers(self, term: str) -> List[str]:
        """Fetch all accession numbers matching a search term."""
        print(f"   🔍 Querying NCBI: {term[:60]}...")
        try:
            handle = Entrez.esearch(
                db="nucleotide",
                term=term,
                rettype="gb",
                retmode="text",
                retmax=100_000, # Large limit to capture full outbreaks
            )
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]
            print(f"      ✅ Found {len(ids)} unique IDs.")
            return ids
        except Exception as e:
            print(f"      ❌ Search failed: {e}")
            return []

    def _filter_and_write(self, records: list, size_thresh_mb: float, output_dir: str, name_key: str, start_idx: int) -> int:
        valid_recs = []
        for r in records:
            size_mb = len(r.seq) / 1_000_000
            if size_mb >= size_thresh_mb:
                acc_id = r.id.split('.')[0] if '.' in r.id else r.id
                r.id = f"{acc_id}_{name_key}"
                r.description = ""
                valid_recs.append(r)

        # Write EACH sequence to a SEPARATE file to ensure it's treated as a distinct strain
        for i, rec in enumerate(valid_recs):
            filepath = Path(output_dir) / f"{name_key}_{start_idx + i}.fasta"
            with open(filepath, "w") as f:
                SeqIO.write([rec], f, "fasta")

        return len(valid_recs)

    def fetch_sequences_chunk(self, id_list: List[str]) -> list:
        """Downloads a batch of sequences with retry logic."""
        attempt = 0
        while attempt < self.max_retries:
            try:
                with Entrez.efetch(
                    db="nucleotide",
                    id=id_list,
                    rettype="fasta",
                    retmode="text",
                ) as handle:
                    return list(SeqIO.parse(handle, "fasta"))
            except KeyboardInterrupt:
                attempt += 1
                if attempt < self.max_retries:
                    print(f"      ⚠️ Download interrupted, retrying ({attempt}/{self.max_retries})...")
                    time.sleep(2 * attempt)
            except Exception as e:
                attempt += 1
                time.sleep(2 * attempt) # Exponential backoff
        return []

    def _auto_detect_genome_size(self, base_query: str) -> float:
        """Fetches up to 20 records to find the median genome size and returns 50% of it."""
        for attempt in range(2):
            try:
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=base_query,
                    rettype="gb",
                    retmode="text",
                    retmax=20,
                )
                record = Entrez.read(handle)
                handle.close()
                ids = record["IdList"]
                if not ids:
                    return 0.0

                with Entrez.efetch(
                    db="nucleotide",
                    id=ids,
                    rettype="fasta",
                    retmode="text",
                ) as handle:
                    recs = list(SeqIO.parse(handle, "fasta"))
                    if not recs: return 0.0
                    lengths = sorted([len(r.seq) for r in recs])
                    median_len = lengths[len(lengths) // 2]
                    return median_len / 1_000_000 # Return median in Mb
            except KeyboardInterrupt:
                print(f"      ⚠️ Auto-detect interrupted (attempt {attempt+1}/2). Retrying...")
                time.sleep(1)
            except Exception as e:
                print(f"      ❌ Auto-detect failed: {e}")
                return 0.0
        print(f"      ⚠️ Auto-detect failed after retries. Proceeding without size filter.")
        return 0.0

    def fetch_and_save_all(
        self,
        query_dict: Dict[str, Tuple[str, float, int]],
        output_folder: str
    ):
        """
        Main Execution Pipeline.
        Args:
            query_dict: {Filename: (Query, SizeMB, MaxCount)}
            output_folder: Path to save.
        """
        out_path = Path(output_folder)
        out_path.mkdir(parents=True, exist_ok=True)
        print(f"📂 Output Directory: {out_path}")

        for name_key, val in query_dict.items():
            search_term, size_thresh, max_count = val[0], float(val[1]), int(val[2])
            target_type = val[3] if len(val) > 3 else "genome"
            if "gene" in search_term.lower():
                target_type = "gene"

            print(f"\n🚀 Processing Task: {name_key}")
            print(f"   📏 Size Threshold: >= {size_thresh} Mb")
            print(f"   🎯 Type: {target_type}")

            # 0. Construct NCBI Query with Native Size Filter & Title Constraints
            # Only apply [title] formatting if it's a simple query (doesn't already contain '[')
            base_query = search_term
            if "[" not in base_query and "]" not in base_query:
                if target_type == "gene":
                    base_query = f'{search_term}'
                else:
                    base_query = f'{search_term}[Organism] AND "complete genome"[title]'

            actual_query = base_query

            # --- AUTO-DETECT GENOME SIZE ---
            max_cutoff_bp = 999999999
            if target_type == "genome" and size_thresh <= 0.0:
                print(f"   🤖 Auto-detecting genome size for: {base_query}...")
                detected_median = self._auto_detect_genome_size(base_query)
                if detected_median > 0:
                    size_thresh = detected_median * 0.7
                    max_cutoff_bp = int(detected_median * 1.5 * 1_000_000)
                    print(f"      ✅ Detected Median Size: {detected_median:.2f} Mb. Setting Safe Bounds: {size_thresh:.2f} - {detected_median * 1.5:.2f} Mb")
                else:
                    print(f"      ⚠️ Auto-detect failed. Proceeding without size filter.")

            if size_thresh and size_thresh > 0:
                cutoff_bp = int(size_thresh * 1_000_000)
                actual_query = f"({base_query}) AND {cutoff_bp}:{max_cutoff_bp}[SLEN]"

            print(f"   🧬 Optimized NCBI Query: {actual_query}")

            # 1. Get IDs (NCBI already filters out small genomes thanks to SLEN)
            acc_ids = self.fetch_accession_numbers(actual_query)
            if not acc_ids:
                print(f"   ⚠️ No results found. Skipping.")
                continue

            # --- OPTIMIZATION: RANDOM SAMPLING ---
            # If user wants 50 genomes but found 5000, we shuffle and take random 50.
            if max_count > 0 and len(acc_ids) > max_count:
                print(f"   📉 Optimization: Sampling {max_count} random genomes from {len(acc_ids)} available.")
                # Randomly select unique IDs without replacement
                acc_ids = random.sample(acc_ids, max_count)
            # -------------------------------------

            # 2. Setup Output File (Clear previous if exists)
            # Find and remove any existing files for this name_key
            for existing_file in out_path.glob(f"{name_key}_*.fasta"):
                existing_file.unlink()

            # 3. Process in Chunks
            total_saved = 0
            chunks = [acc_ids[i:i + self.chunk_size] for i in range(0, len(acc_ids), self.chunk_size)]

            with tqdm(total=len(chunks), desc="   ⬇️ Downloading", unit="chunk") as pbar:
                for chunk in chunks:
                    raw_seqs = self.fetch_sequences_chunk(chunk)
                    if raw_seqs:
                        count = self._filter_and_write(raw_seqs, size_thresh, str(out_path), name_key, total_saved)
                        total_saved += count
                    pbar.update(1)

            if total_saved > 0:
                print(f"   ✅ Success: {total_saved} genomes saved (separated for strain-aware processing)")
            else:
                print(f"   ⚠️ Warning: 0 genomes met the size criteria.")

        print("\n✨ All download tasks completed.")
