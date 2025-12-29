import os
import ssl
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from Bio import Entrez, SeqIO
from tqdm import tqdm

class SequenceFetcher:
    """
    A robust NCBI downloader that handles large datasets, network retries,
    and dynamic genome size filtering.
    """
    def __init__(
        self, 
        email: str, 
        api_key: Optional[str] = None, 
        chunk_size: int = 200, 
        max_retries: int = 3
    ):
        """
        Args:
            email: User email (required by NCBI).
            api_key: NCBI API Key (optional, boosts speed to 10 reqs/sec).
            chunk_size: Number of IDs to fetch per batch (default 200).
            max_retries: Number of retry attempts for failed chunks.
        """
        self.email = email
        self.chunk_size = chunk_size
        self.max_retries = max_retries
        
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            
        # Bypass SSL verification issues on some legacy systems
        self.ssl_context = ssl._create_unverified_context()

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
                ssl_context=self.ssl_context
            )
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]
            print(f"      ✅ Found {len(ids)} unique IDs.")
            return ids
        except Exception as e:
            print(f"      ❌ Search failed: {e}")
            return []

    def _filter_and_write(self, records: list, size_thresh_mb: float, filepath: Path) -> int:
        """
        Filters records by size and appends valid ones to the file immediately.
        Returns the count of valid sequences written.
        """
        if size_thresh_mb is None:
            valid_recs = records
        else:
            cutoff_bp = size_thresh_mb * 1_000_000
            valid_recs = [r for r in records if len(r.seq) >= cutoff_bp]

        if not valid_recs:
            return 0

        # Append to file (Incremental Write)
        with open(filepath, "a") as f:
            SeqIO.write(valid_recs, f, "fasta")
            
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
                    ssl_context=self.ssl_context
                ) as handle:
                    return list(SeqIO.parse(handle, "fasta"))
            except Exception as e:
                attempt += 1
                time.sleep(2 * attempt) # Exponential backoff
        return []

    def fetch_and_save_all(
        self, 
        query_dict: Dict[str, Tuple[str, float]], 
        output_folder: str
    ):
        """
        Main Execution Pipeline.
        
        Args:
            query_dict: Dictionary mapping Filenames -> (SearchQuery, SizeThresholdMB)
            output_folder: Directory to save FASTA files.
        """
        out_path = Path(output_folder)
        out_path.mkdir(parents=True, exist_ok=True)
        print(f"📂 Output Directory: {out_path}")
        
        for name_key, (search_term, size_thresh) in query_dict.items():
            print(f"\n🚀 Processing Task: {name_key}")
            print(f"   📏 Size Threshold: >= {size_thresh} Mb")
            
            # 1. Get IDs
            acc_ids = self.fetch_accession_numbers(search_term)
            if not acc_ids:
                print(f"   ⚠️ No results found. Skipping.")
                continue

            # 2. Setup Output File (Clear previous if exists)
            filename = f"{name_key}.fasta"
            file_path = out_path / filename
            if file_path.exists():
                file_path.unlink() # Start fresh

            # 3. Process in Chunks (Memory Safe)
            total_saved = 0
            chunks = [acc_ids[i:i + self.chunk_size] for i in range(0, len(acc_ids), self.chunk_size)]
            
            # TQDM Progress Bar
            with tqdm(total=len(chunks), desc="   ⬇️ Downloading", unit="chunk") as pbar:
                for chunk in chunks:
                    # Download
                    raw_seqs = self.fetch_sequences_chunk(chunk)
                    
                    if raw_seqs:
                        # Filter & Write immediately to disk
                        count = self._filter_and_write(raw_seqs, size_thresh, file_path)
                        total_saved += count
                    
                    pbar.update(1)

            # 4. Summary
            if total_saved > 0:
                print(f"   ✅ Success: {total_saved} genomes saved to {filename}")
            else:
                print(f"   ⚠️ Warning: 0 genomes met the size criteria.")

        print("\n✨ All download tasks completed.")