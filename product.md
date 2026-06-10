# Rational Primer Design System — Product & Architecture Document

## 1. Overview

**Rational Primer Design** is a bioinformatics pipeline for *in-silico* design, validation, and multiplex assembly of PCR/qPCR primers and TaqMan probes at population scale. It uses a K-mer mining algorithm over whole genomes, in-silico PCR simulation, and optional LLM-based cognitive evaluation.

**Purpose:** Given a target species and a panel of background (exclusivity) species, produce a set of primer pairs and probes that amplify the target with high sensitivity and do not cross-react with the background.

**Version:** 1.0.3  
**Language:** Python 3.11+ (backend), TypeScript 6 / Next.js 16 (frontend)  
**License:** MIT  

---

## 2. System Architecture (Three-Tier)

```
┌─────────────────────────────────────────────────────────────┐
│                    CLIENT LAYER                              │
│  ┌───────────┐  ┌───────────┐  ┌───────────┐              │
│  │ Web UI    │  │ CLI       │  │ REST API  │              │
│  │ (Next.js) │  │ (argparse)│  │ (curl/http)│             │
│  └─────┬─────┘  └─────┬─────┘  └─────┬─────┘              │
└────────┼───────────────┼───────────────┼────────────────────┘
         │               │               │
         ▼               ▼               ▼
┌─────────────────────────────────────────────────────────────┐
│                    API LAYER (FastAPI)                        │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐                  │
│  │ Job Mgmt │  │ AI Chat  │  │ File     │                  │
│  │ /api/jobs│  │ /api/ai  │  │ Browser  │                  │
│  └─────┬────┘  └────┬─────┘  └──────────┘                  │
│        │            │                                       │
│  ┌─────┴────────────┴──────────────────────────────────┐   │
│  │ subprocess.Popen  │  PROCESS_REGISTRY dict           │   │
│  └──────────────────────────────────────────────────────┘   │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│                   PIPELINE LAYER (Python)                     │
│                                                               │
│  ┌──────────┐   ┌────────────┐   ┌──────────────┐          │
│  │ Fetcher  │──▶│ Constructor│──▶│  Designer    │          │
│  │ (NCBI)   │   │ (Library)  │   │  (K-mer)     │          │
│  └──────────┘   └────────────┘   └──────┬───────┘          │
│                                         │                    │
│  ┌──────────┐   ┌────────────┐   ┌──────▼───────┐          │
│  │ Prober   │◀──│ Validator  │◀──│ Candidate    │          │
│  │ (TaqMan) │   │ (In-silico │   │ Pairs        │          │
│  └────┬─────┘   │  PCR)      │   └──────────────┘          │
│       │         └────────────┘                              │
│       │         ┌────────────┐                              │
│       └────────▶│   Final    │                              │
│                 │   Assay    │                              │
│                 │   + BLAST  │                              │
│                 └────────────┘                              │
│                                                               │
│  ┌──────────┐   ┌────────────┐   ┌──────────────┐          │
│  │  AI Core │   │ Multiplex  │   │  Advanced    │          │
│  │ (LLM)    │   │ Assembler  │   │  PCR Engine  │          │
│  └──────────┘   └────────────┘   └──────────────┘          │
└─────────────────────────────────────────────────────────────┘
```

---

## 3. Pipeline Flow (Six Stages)

### Stage 0 — NCBI Data Fetching
**Module:** `rational_design/fetcher.py` — `SequenceFetcher` class

```
User query (species name, type=genome/gene, count, size threshold)
    │
    ▼
NCBI ESearch (term + SLEN filter after genome size auto-detect)
    │
    ▼
Random sampling (if count < available IDs)
    │
    ▼
Chunked EFetch (10 at a time, retry with exponential backoff)
    │
    ▼
Size filter → write per-strain FASTA files to 0_raw_data/
```

- Genome size auto-detect: median of 20 random records
- SLEN filter: `lower = median * 0.7`, `upper = median * 1.5`
- Chunk size: 10, sleep 3.4s between chunks (NCBI rate limit)

### Stage 1 — Dataset Construction
**Module:** `rational_design/constructor.py` — `LibraryConstructor` class

```
Raw FASTA files (0_raw_data/target/, 0_raw_data/background/)
    │
    ▼
Group by filename prefix (e.g., t1, b1, b2) or treat as single pool
    │
    ▼
Random sampling per group (design_target_sampling_size, etc.)
    │
    ▼
Split into design set and validation set (separate files)
    │
    ▼
Output: 1_workspace/design/{target,background}.fasta
        1_workspace/validate/{target,background}.fasta
```

- Record ID convention: `strain|contig` (pipe-separated)
- 1 FASTA file = 1 biological strain, regardless of contig count

### Stage 2 — Primer Design (K-mer Mining)
**Module:** `rational_design/designer.py` — `PrimerDesigner` class

This is the core algorithm. It runs in three phases per K-mer length (18–22 bp):

```
Target genomes (FASTA)                    Background genomes (FASTA)
    │                                            │
    ▼                                            ▼
┌─────────────────────────────────────────────────────────────┐
│  PHASE 1: Exhaustive K-mer Mining on Target                  │
│                                                               │
│  For each target strain (multiprocessing.Pool):               │
│    stride=1 sliding window over entire genome                 │
│    2-bit bit-packed integer hashing (A=00, C=01, G=10, T=11) │
│    Filter each window:                                        │
│      ├─ is_bio_safe():   GC content 40-75%, no homopolymer≥5 │
│      ├─ calculate_entropy():  Shannon entropy ≥ 1.4           │
│      └─ is_low_complexity_sequence():  no tandem repeats      │
│                                                               │
│  Cross-strain conservation counting (Counter):                 │
│    Keep k-mers present in ≥ design_min_conservation fraction  │
│    Optionally expand below-threshold k-mers via IUPAC variants│
│    Cap at 50,000 highest-conservation k-mers → "elite set"   │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  PHASE 2: Background Filtering                                │
│                                                               │
│  Scan each background genome for elite k-mer presence:        │
│    stride=1, return set of found hashes per strain            │
│    Counter accumulates strain-level counts across background  │
│                                                               │
│  Mode selection (config.background_filter_mode):              │
│    "strict"       → keep only k-mers with zero background hits│
│    "amplicon_level" → always use rescue scoring               │
│    "auto" (default) → strict mode, but if clean pool < 500   │
│                        or clean ratio < 0.05 → rescue mode    │
│                                                               │
│  Rescue scoring (per k-mer):                                  │
│    score = 100 × target_freq - 300 × background_freq          │
│    Filter: drop k-mers with background_freq > 0.25            │
│    Keep top 5000 scored k-mers for pairing                    │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  PHASE 3: Primer Pair Formation                               │
│                                                               │
│  Reference genome (first target strain):                       │
│    Rolling hash to locate anchor positions of all elite k-mers│
│    For each anchor at position A:                              │
│      Look forward for second anchor at position B              │
│      product_size = B - A + k  ∈ [min_p, max_p]               │
│      Forward primer = seq[A..A+k]                              │
│      Reverse primer = rev_comp(seq[B..B+k])                    │
│      Check Tm (55-68°C), Tm delta (≤3°C)                      │
│      Check structural QC (hairpin, dimer)                      │
│      Keep up to design_max_candidates pairs per K length       │
│                                                               │
│  Rescue mode extra step:                                       │
│    For each candidate pair, run in-silico PCR on background    │
│    genomes using find_primer_hits(). Reject if valid amplicon  │
│    is produced in background (final_max_background_amplicon=0) │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  FINAL: Scoring & Export                                      │
│                                                               │
│  Multithreaded prevalence scoring across all target strains:   │
│    prevalence = % strains with ≥1 hit                          │
│    avg_copy = mean copy number across positive strains         │
│  Sort by avg copy number, keep top design_max_candidates      │
│  Export to 2_candidates/candidates.csv                        │
└─────────────────────────────────────────────────────────────┘
```

**Key data structures:**

| Structure | Type | Description |
|-----------|------|-------------|
| K-mer hash | `int` (64-bit) | 2-bit encoded, max 22 bp → 44 bits |
| `target_pool` | `Counter[int]` | Strain-level conservation counts |
| `valid_kmers` | `dict[int, int]` | Conservation-filtered k-mers (hash → count) |
| `elite_targets` | `dict[int, int]` | Top 50K k-mers with counts |
| `bg_hit_counter` | `Counter[int]` | Background strain counts per k-mer |
| Anchor set | `set[int]` | K-mers passed for pairing (strict or rescue) |
| Candidate pairs | `list[dict]` | `{Fwd: str, Rev: str, Pos: int}` |

### Stage 3 — In-Silico Validation
**Module:** `rational_design/validator.py` — `InSilicoValidator` class

```
Primer candidates (FASTA)      Target/Background genomes (FASTA)
    │                                    │
    ▼                                    ▼
┌─────────────────────────────────────────────────────────────┐
│  Pigeonhole Segment Search Algorithm:                        │
│                                                               │
│  For non-degenerate primers:                                  │
│    Split primer into (max_mismatch+1) segments                │
│    At least one segment must match exactly (pigeonhole)       │
│    Use Python str.find() for fast exact segment matching      │
│    For each candidate start position, compute Hamming dist    │
│                                                               │
│  For degenerate primers (IUPAC):                              │
│    Check every position in genome (no pigeonhole acceleration)│
│    Use _iupac_mismatch_count() for weighted mismatch counting │
│                                                               │
│  Orientation: check both Fwd + Rev_rc and Rev + Fwd_rc        │
│  Product: forward.pos < reverse.pos, size ∈ [min, max]        │
│                                                               │
│  Multiprocessing: batched by strain group (batch_size varies  │
│    by available RAM: 5-100 strains per batch)                 │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
Output files:
  results_target.csv       — per-strain PCR products for target
  results_background.csv   — per-strain PCR products for background
  pcr_results_summary.csv  — sensitivity %, copy number, consensus
```

**Validation outputs per candidate pair:**

| Column | Description |
|--------|-------------|
| `Strain` | Strain ID (filename stem or contig prefix) |
| `Primer Pair` | `Set_1`, `Set_2`, etc. |
| `PCR Product Sequence` | Amplicon sequence (if positive) |
| `Total_Bands` | Total amplicons detected in strain |
| `Mismatch_Map` | `F:[map]; R:[map]` showing mismatch positions |
| `Forward Primer` / `Reverse Primer` | Primer sequences |
| `F_Bind_Seq` / `R_Bind_Seq` | Actual binding site sequences |

**Auto-retry loop (cli.py):** Three cycles with progressively relaxed constraints:
1. Original parameters
2. Conservation -0.05, sensitivity -5%
3. Conservation -0.10, sensitivity -10%

### Stage 4 — Probe Design
**Module:** `rational_design/prober.py` — `ProbeSelector` class

```
Input: pcr_results_summary.csv (at sensitivity threshold)
       results_target.csv (amplicon sequences)
       results_background.csv (optional)
    │
    ▼
For each primer pair:
  Collect all unique amplicon sequences from target strains
  Compute specificity = 1 - (background_hits / total_background)
    │
    ▼
Probe search (per amplicon):
  Search zone: exclude 15bp from each primer end
  Sliding window 20-25bp
  Filter:
    ├─ Tm = primer_Tm + 5~12°C (SantaLucia 1998, DNA_NN3)
    ├─ GC content 35-65%
    ├─ No starting G at 5' end
    ├─ No homopolymer runs ≥4
    └─ Must bind in ALL amplicon variants
    │
    ▼
Cross-check with primers:
  No hairpin ≥4bp
  No 3' self-dimer ≥4bp
  No 3' cross-dimer ≥4bp with either primer
    │
    ▼
Score: closest to (primer_Tm + 8°C) and GC = 50%
Select highest-scoring probe
    │
    ▼
Output: FINAL_ASSAY.csv
  Set_ID, Forward Primer, Reverse Primer, Probe_Seq, Probe_Tm,
  Amplicon_Size, Sensitivity, Specificity, Target_Gene
```

### Stage 5 — BLAST Annotation (Optional)
**Module:** `rational_design/prober.py` — `ProbeSelector.run_blast_annotation()`

- Runs `blastx` against NCBI NR database (30s timeout per query)
- Extracts `Target_Gene` name from top hit description
- Annotates `FINAL_ASSAY.csv` with `Target_Gene` column

### Stage 6 — Reporting & Audit Trail
**Module:** `rational_design/cli.py` — `build_validation_report()`, `write_audit_trail()`

- Generates `4_validation_report/` with formatted assay details
- Cross-contamination traceback report
- Immutable `audit_trail.json` with timestamp, parameters, results, AI decisions

---

## 4. Module Architecture (Detailed)

### 4.1 Core Pipeline Modules

| Module | Class | Entry Point | Key Dependencies |
|--------|-------|-------------|------------------|
| `fetcher.py` | `SequenceFetcher` | `fetch_and_save_all()` | Bio.Entrez, requests |
| `constructor.py` | `LibraryConstructor` | `construct()` | Bio.SeqIO |
| `designer.py` | `PrimerDesigner` | `design()` | multiprocessing, Bio.Seq, Bio.SeqUtils |
| `validator.py` | `InSilicoValidator` | `validate()` | Levenshtein, Bio.SeqIO, multiprocessing |
| `prober.py` | `ProbeSelector` | `design()` | Bio.SeqUtils, pandas, requests (BLAST) |

### 4.2 Support Modules

| Module | Class/Function | Purpose |
|--------|----------------|---------|
| `multiplex.py` | `check_dimer()`, `check_hairpin()`, `evaluate_multiplex_kit()` | Cross-primer compatibility, multiplex assembly |
| `insilico_pcr_advanced.py` | `IndustrialEngine` | High-performance PCR for known/custom primers |
| `ai_core.py` | `AIExpertAgent`, `AssayEvaluator` | LLM integration, cognitive evaluation |
| `utils.py` | `DualLogger`, `generate_iupac_consensus()` | I/O, Tm, GC computation |
| `analytics.py` | `ResultAnalysisEngine` | Historical run indexing and analysis |

### 4.3 Interface Modules

| Module | Class/Function | Purpose |
|--------|----------------|---------|
| `cli.py` | `run_full_pipeline()` | Command-line pipeline orchestrator |
| `term.py` | Wizard functions | Interactive Rich terminal interface |
| `api.py` | FastAPI app factory | REST API server |
| `api_server.py` | `main()` | uvicorn launcher |
| `web_jobs.py` | `run_*()` functions | Background job runners for API |

### 4.4 Module Dependency Graph

```
cli.py / api.py
  ├── fetcher.py
  ├── constructor.py
  ├── designer.py
  │    ├── utils.py (via IUPAC consensus)
  │    └── validator.py (find_primer_hits import)
  ├── validator.py
  │    ├── utils.py (nuclear_ram_flush)
  │    └── Levenshtein (external)
  ├── prober.py
  │    ├── utils.py (IUPAC consensus)
  │    ├── multiplex.py (dimer, hairpin checks)
  │    └── validator.py (IUPAC table)
  ├── multiplex.py
  ├── ai_core.py
  │    ├── utils.py (generate_batch_analytical_summary)
  │    └── openai (external)
  └── analytics.py
```

---

## 5. Configuration System

### Parameter Hierarchy

```
Hard-coded defaults (cli.py lines 165-185)
    │
    ▼
config/parameters.json (user-editable JSON file)
    │
    ▼
CLI arguments (--params, --out, --email, etc.)
    │
    ▼
AI proposal overrides (via JSON proposal block)
```

### Parameter Categories

**Design set:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `design_target_sampling_size` | 50 | Max target genomes for K-mer mining (0=all) |
| `design_background_sampling_size` | 50 | Max background genomes for K-mer blacklisting |
| `primer_length_min` | 18 | Minimum K-mer/primer length (bp) |
| `primer_length_max` | 22 | Maximum K-mer/primer length (bp) |
| `primer_tm_min` | 55.0 | Minimum primer melting temperature (°C) |
| `primer_tm_max` | 68.0 | Maximum primer melting temperature (°C) |
| `design_min_conservation` | 0.75 | Min fraction of target strains containing a K-mer |
| `design_max_candidates` | 50 | Max candidate pairs to output per K-length |

**Validation set:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `validation_target_sampling_size` | 0 | Max target genomes for validation (0=all) |
| `validation_background_sampling_size` | 0 | Max background genomes for validation |
| `min_sensitivity` | 90.0 | Minimum acceptable sensitivity (%) |
| `max_mismatch` | 2 | Max mismatches for in-silico PCR tolerance |

**General:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `product_size_min` | 120 | Minimum amplicon size (bp) |
| `product_size_max` | 400 | Maximum amplicon size (bp) |
| `cpu_cores` | 0 | CPU cores (0 = auto-detect all) |
| `degenerate_primers` | true | Enable IUPAC degenerate primer generation |
| `max_iupac_per_primer` | 2 | Max degenerate positions per primer |
| `enable_blast` | true | Enable BLASTX gene annotation |
| `auto_relax_constraints` | true | Auto-relax thresholds if first pass finds 0 candidates |

**Adaptive background filtering (new in v1.0.3):**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `background_filter_mode` | "auto" | `strict`, `amplicon_level`, or `auto` |
| `strict_min_clean_kmers` | 500 | Min clean k-mers to stay in strict mode |
| `strict_min_clean_ratio` | 0.05 | Min clean/elite ratio to stay in strict mode |
| `rescue_top_kmers_for_pairing` | 5000 | Max k-mers to keep in rescue mode |
| `rescue_max_single_primer_background_freq` | 0.25 | Max background freq for a rescue k-mer |
| `rescue_background_penalty_weight` | 300.0 | Penalty multiplier for background hits in scoring |
| `final_max_background_amplicon_hits` | 0 | Max allowed background amplicons per pair |

### Config Files on Disk

```
config/
├── parameters.json         # Core pipeline parameters
├── background.json         # NCBI background queries (key → [term, size, count, type])
├── target.json             # NCBI target queries (same format)
├── b_conf_shared.json      # Shared background for multiplex runs
├── t_conf.json             # Multiplex target config (single target)
├── t_conf_target1-4.json   # Per-target configs for multiplex
├── user_settings.json      # Saved user settings (email, etc.)
└── ai_settings.json        # Runtime AI endpoint and model name
```

---

## 6. Key Algorithms (Detailed)

### 6.1 Bit-Packed K-mer Hashing

DNA bases are encoded as 2-bit integers:
- A = `00`, C = `01`, G = `10`, T = `11`
- A 22-mer fits in 44 bits (fits in a Python `int`)
- Rolling hash: `h = ((h << 2) | val) & mask` for O(1) sliding window
- Pre-computed `DNA_MAP` bytearray for O(1) char-to-code lookup

### 6.2 Pigeonhole Segment Search (Validator)

For a primer of length *n* with tolerance *m* mismatches:
1. Split primer into *m+1* segments (as evenly as possible)
2. By the pigeonhole principle, at least one segment must match exactly
3. Use Python's fast C-level `str.find()` for exact segment matching
4. For each candidate start, compute Hamming distance (or IUPAC-weighted mismatch)

This gives ~10-100× speedup over brute-force sliding window.

### 6.3 Tm Calculation (SantaLucia 1998)

- Nearest-neighbor thermodynamics with `DNA_NN3` table
- Salt correction: Na=50mM, Mg=3.0mM, dNTPs=0.8mM
- Oligo concentration: 250nM
- For degenerate primers: Tm approximated from non-IUPAC bases only

### 6.4 Multiplex Kit Scoring

Score formula (0-100):
```
score = 100
  - tm_penalty          (Tm span > 4°C)
  - total_dimer_score   (cross-dimer weighted sum)
  - hairpin_penalty     (any hairpin ≥ 4bp)
  - self_dimer_penalty  (any self-dimer ≥ 4bp 3')
  - probe_penalty       (missing probe in qPCR mode)
  - size_penalty        (amplicon overlap in gel mode)
```

Verdicts: EXCELLENT (≥85), GOOD (70-84), MARGINAL (50-69), POOR (<50)

---

## 7. Multiplex Assembly Architecture

### Two Modes

**1. Sequential Design Mode** (for new multiplex kits):
- Each target designed independently through full pipeline
- Other targets injected as additional background strains
- `shared_primer_context.json` accumulates accepted assays
- Gate 4 check (Python + AI): cross-target Tm span ≤4°C, no 3' cross-dimers ≥4bp

**2. Combinatorial Analysis Mode** (for existing assays):
- Load existing `FINAL_ASSAY.csv` files from individual target folders
- Extract top 5 candidates per target (avoids combinatorial explosion)
- Cartesian product → scored via `evaluate_multiplex_kit()`
- Return top K kits sorted by score

---

## 8. AI Integration Architecture

### Two AI Roles

**1. AssayEvaluator** (batch-level gate keeper):
```
Input: Pre-computed analytical summary (Tm, GC, QC flags)
  │
  ▼
Gate 1: Tm range OK? Tm Delta ≤3°C?
Gate 2: QC flags all pass?
Gate 3: Sensitivity ≥ min_sensitivity, GC 40-60%?
Gate 4: (multiplex only) Cross-target Tm span ≤4°C?
  │
  ▼
Output: {overall_verdict, next_action, confidence_score, top_3_assays}
  Actions: ACCEPT_AND_STOP | REJECT_AND_CONTINUE
```

- Called between validation batches (every 50 candidates)
- Sub-batch retry: if rejected, try next 10 candidates (up to 3 times)

**2. AIExpertAgent** (chatbot / configuration proposer):
- Available via web API, terminal, or resident sidebar
- System prompts include: molecular biology rules, hardware-aware sampling limits, 4 proposal JSON schemas
- Can propose: `propose_design`, `propose_local_design`, `propose_validation`, `propose_multiplex`
- Auto-run capability: when `run_immediately: true`, API creates and starts job directly
- Proposal extraction: regex parsing for ```json ... ``` fences

**LLM Backend:** OpenAI-compatible (`/v1/chat/completions` with JSON mode)
- Default: Ollama at `http://localhost:11434/v1`
- Model auto-detection via Ollama `/api/tags`

---

## 9. Frontend Architecture

### Stack
- **Framework:** Next.js 16.2.6 (App Router, static generation)
- **UI Library:** React 19.2.6
- **Language:** TypeScript 6.0.3

### Routes (all served by single Dashboard component)

| Route | View (`view` prop) | Purpose |
|-------|-------------------|---------|
| `/` | — | Landing page |
| `/dashboard` | `dashboard` | Job history, results, live monitoring |
| `/design/local` | `local` | Local FASTA file-based design |
| `/design/auto` | `auto` | NCBI keyword-based design |
| `/ai` | `ai` | AI chat + proposal interface |
| `/validate` | `validate` | Known primer validation |
| `/multiplex` | `multiplex` | Multiplex design & analysis |
| `/about` | `about` | Static information |

### Component Structure

```
Dashboard (root component, ~1600 lines in Dashboard.tsx)
  ├── Header (brand, navigation, status)
  ├── SharedConfigPanel (language, NCBI email, AI endpoint, params)
  │    ├── Parameter inputs (number/boolean)
  │    └── System status display
  ├── RecentJobs (recent job list with status badges)
  └── Workspace Area (view-dependent)
       ├── DashboardHome
       │    ├── JobList
       │    ├── ResultsPanel
       │    └── MonitorPanel (terminal emulator, logs, cancel/delete)
       ├── AutoMode
       │    ├── QueryEditor (target keywords)
       │    ├── QueryEditor (background keywords)
       │    ├── Summary metrics
       │    └── Run button
       ├── LocalMode
       │    ├── File path inputs
       │    └── Upload toggle
       ├── AiMode
       │    ├── ChatMessages
       │    ├── ProposalSummary (structured proposal display)
       │    └── Run proposal button
       ├── ValidationMode
       │    ├── Primer pair inputs
       │    └── Mode toggle (local file / NCBI auto)
       └── MultiplexModePanel
            ├── Target folder selection
            ├── Assay type (qPCR/gel)
            └── Analyze existing button
```

### State Management
- **Global state:** `lang`, `parameters`, `jobs`, `job`, `logs`, `results`, `status`
- **Parameter defaults:** Fetched from `GET /api/default-parameters` on mount
- **Persistence:** `localStorage` for language, NCBI email, AI settings, file paths

### API Client
```typescript
async function api<T>(url: string, options?: RequestInit): Promise<T>
```
- Single generic wrapper around `fetch()`
- JSON body serialization
- Error response parsing with detail extraction
- Used for all API interactions

### Localization
- 300+ i18n keys per language (English and Vietnamese)
- Keys organized by UI section: `nav`, `dashboard`, `sidebar`, `design`, `validate`, `about`, `common`

---

## 10. REST API Design

### Base URL: `http://host:port/api`

### Endpoint Summary

| Method | Path | Purpose | Request Body | Response |
|--------|------|---------|-------------|----------|
| GET | `/health` | Health check | — | `{status: "ok"}` |
| GET | `/status` | System status | — | `{disk, jobs, ai_active, models[]}` |
| GET | `/default-parameters` | Default pipeline params | — | `dict` of parameters |
| GET | `/ai/models` | Available LLM models | — | `{models: [{id, name}]}` |
| GET | `/files/browse` | File browser | Query: `path`, `kind` | `{current, parent, entries[]}` |
| GET | `/jobs` | List all web jobs | Query: `page`, `limit` | `[{id, status, created_at, ...}]` |
| GET | `/history` | Legacy + web run history | — | `[run objects]` |
| POST | `/jobs/local` | Start local design | `LocalPipelineRequest` | `{id, status}` |
| POST | `/jobs/auto` | Start NCBI auto design | `AutoPipelineRequest` | `{id, status}` |
| POST | `/jobs/validate` | Validate known primers | `ValidationRequest` | `{id, status}` |
| POST | `/jobs/multiplex/local` | Local multiplex | `LocalMultiplexRequest` | `{id, status}` |
| POST | `/jobs/multiplex/auto` | NCBI auto multiplex | `AutoMultiplexRequest` | `{id, status}` |
| POST | `/jobs/multiplex/analyze` | Analyze existing | `MultiplexAnalyzeRequest` | `{id, status}` |
| POST | `/jobs/upload` | Upload + run | Multipart form | `{id, status}` |
| GET | `/jobs/{id}` | Job status | — | `{id, status, ...}` |
| GET | `/jobs/{id}/logs` | Job logs | — | text/plain |
| GET | `/jobs/{id}/results` | Job output files | — | `{files: [{name, path, size}]}` |
| POST | `/jobs/{id}/cancel` | Cancel job | — | `{status: "cancelled"}` |
| GET | `/jobs/{id}/archive` | ZIP download | — | application/zip |
| DELETE | `/jobs/{id}` | Delete job | — | `{status: "deleted"}` |
| POST | `/ai/chat` | AI chatbot | `ChatRequest` | `{reply, proposal?, error?}` |

### Job Lifecycle
```
pending → running → completed
                → failed
                → cancelled
```

### Background Job Execution
1. API receives request → writes `request.json` to job directory
2. Spawns `subprocess.Popen(["python", "-m", "rational_design.web_jobs", command, "--request", request_path])`
3. Stores `Popen` object in `PROCESS_REGISTRY` dict (keyed by job ID)
4. Client polls `GET /api/jobs/{id}` for status updates
5. Status inferred from output files: `FINAL_ASSAY.csv` → completed, `error.log` → failed

---

## 11. Data Flow Diagram (End-to-End)

```
User Input (species keywords / local FASTA paths)
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  CONFIGURATION                                               │
│  parameters.json + AI proposal + user overrides              │
│  → merged dict passed through pipeline                       │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  DATA ACQUISITION                                            │
│  SequenceFetcher: NCBI ESearch → EFetch → per-strain FASTA   │
│  LibraryConstructor: sample → split → strain|contig FASTA    │
│                                                               │
│  Output files:                                                │
│    1_workspace/design/target.fasta                            │
│    1_workspace/design/background.fasta                        │
│    1_workspace/validate/target.fasta                          │
│    1_workspace/validate/background.fasta                      │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  PRIMER DESIGN (loop k=18..22)                               │
│  PrimerDesigner.design()                                      │
│                                                               │
│  Input:  design/target.fasta, design/background.fasta         │
│  Internal:  target_strains: dict[str, list[str]]              │
│             bg_strains: dict[str, list[str]]                   │
│  Algorithm:                                                    │
│    Phase 1 (K-mer mining):                                    │
│      target_strains → Counter(hash) → conservation filter    │
│      → elite_targets: dict[hash, count]                      │
│    Phase 2 (Background screen):                               │
│      elite_targets + bg_strains → bg_hit_counter: Counter    │
│      → select_anchors_for_pairing() → anchors: set[hash]     │
│    Phase 3 (Pair formation):                                  │
│      anchors + reference genome → candidate_pairs: list[dict] │
│      rescue mode: filter_pairs_by_background_amplicon()       │
│    Scoring:                                                   │
│      candidate_pairs + all_target_strains → scored CSV        │
│                                                               │
│  Output: 2_candidates/candidates.csv                          │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  IN-SILICO VALIDATION (iterative batches)                     │
│  InSilicoValidator.validate()                                  │
│                                                               │
│  Input: validate/target.fasta, validate/background.fasta,     │
│         candidates.fasta (converted from candidates.csv)      │
│  Algorithm: Pigeonhole segment search + Hamming distance      │
│  Multiprocessing: batched by strain (RAM-aware batch size)    │
│                                                               │
│  Output:                                                      │
│    results_target.csv     — per-strain PCR hits               │
│    results_background.csv — per-strain background hits        │
│    pcr_results_summary.csv — sensitivity% per candidate       │
└─────────────────────────────────────────────────────────────┘
    │
    ├── AI Evaluation (optional, between batches)
    │   AssayEvaluator.evaluate_candidates() → ACCEPT/REJECT
    │   Sub-batch retry on rejection
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  PROBE DESIGN (if qPCR mode)                                  │
│  ProbeSelector.design()                                       │
│                                                               │
│  Input: pcr_results_summary.csv + results_target.csv          │
│  Algorithm: Sliding window on amplicon, Tm/GC/dimer filters   │
│                                                               │
│  Output: FINAL_ASSAY.csv (with Probe_Seq, Specificity)        │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  BLAST ANNOTATION (optional)                                  │
│  ProbeSelector.run_blast_annotation()                         │
│  blastx on NCBI NR → Target_Gene column in FINAL_ASSAY.csv   │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│  REPORTING & AUDIT                                           │
│  build_validation_report() → 4_validation_report/             │
│  traceback_cross_contamination() → cross_contamination.json   │
│  write_audit_trail() → audit_trail.json (immutable)          │
└─────────────────────────────────────────────────────────────┘
```

---

## 12. Output File Structure (per run)

```
runs/
└── YYYYMMDD-HHMMSS-{uuid}/
    └── primer-design-run/
        ├── 0_raw_data/
        │   ├── target/
        │   │   └── {strain_id}.fasta  (per-strain NCBI downloads)
        │   └── background/
        │       └── {strain_id}.fasta
        ├── 1_workspace/
        │   ├── design/
        │   │   ├── target.fasta       (sampled, merged for design)
        │   │   └── background.fasta   (sampled, merged for design)
        │   └── validate/
        │       ├── target.fasta       (separate set for validation)
        │       └── background.fasta
        ├── 2_candidates/
        │   └── candidates.csv         (raw primer pairs from designer)
        ├── 3_validation/
        │   ├── results_target.csv     (per-strain target PCR results)
        │   ├── results_background.csv (per-strain background PCR)
        │   ├── pcr_results_summary.csv (sensitivity, copy, consensus)
        │   ├── master_results_target.csv (accumulated across batches)
        │   └── master_results_background.csv
        ├── 4_validation_report/
        │   └── (formatted reports)
        ├── FINAL_ASSAY.csv            (accepted assays with probes)
        ├── audit_trail.json           (immutable run record)
        ├── ai_expert_report.json      (AI evaluation report)
        ├── cross_contamination_report.json
        ├── shared_primer_context.json (for multiplex runs)
        └── pipeline_log.txt           (tee'd stdout/stderr)
```

---

## 13. Testing

### Test Files

| File | Tests | Coverage |
|------|-------|----------|
| `tests/test_regressions.py` | 10 | Pigeonhole hit finding, product masking, traceback (3 variants), multiplex rejection, assembly dedup, constructor sampling, low-complexity detection, prober specificity |
| `tests/test_api.py` | 8 | Default params, missing paths, local/auto/validation/multiplex job creation, AI proposal extraction |
| `tests/test_ai_core.py` | 2 | System prompt content verification |

### Test Data
- `test_data/target/1.fasta` — small test genome
- `test_data/background/1.fasta` — small test genome

### Running Tests
```bash
pytest tests/ -v
```

---

## 14. Deployment

### Requirements
- Python 3.9+
- Node.js 20+ (for frontend)
- NCBI Entrez API access (requires email)
- Optional: Ollama or OpenAI-compatible LLM server for AI features

### Installation
```bash
pip install -e .
cd frontend && npm install
```

### Running
```bash
# CLI mode
rational-design --local-target ./target --local-bg ./background --out ./output

# API server
rational-design-api

# Frontend (separate terminal)
cd frontend && npm run dev
```

### Configuration
- Copy `config/parameters.json` and edit parameters
- Set NCBI email in web UI sidebar or `config/user_settings.json`
- For AI features: set AI endpoint in web UI sidebar

---

## 15. System Requirements & Limits

- **RAM:** 16 GB minimum for bacterial genomes (8 GB for very small targets)
- **CPU:** Multi-core recommended (auto-detects all cores)
- **Disk:** Varies with genome count (~1 GB per 100 bacterial genomes)
- **K-mer memory:** 50,000 elite k-mers × 8 bytes ≈ 400 KB (ignoring overhead)
- **Background index:** Scales linearly with background genome count
- **NCBI rate limit:** 10 requests/second, 3.4s sleep between EFetch chunks

### RAM-Based Sampling Limits (AI Auto-Config)

| RAM | Gram+ Bacteria | Gram- Bacteria | Fungi | Virus |
|-----|---------------|---------------|-------|-------|
| 8 GB | 20 | 15 | 10 | 500 |
| 16 GB | 70 | 50 | 20 | 1000 |
| 32 GB | 150 | 100 | 50 | 2000 |

---

## 16. Glossary

| Term | Definition |
|------|------------|
| K-mer | Fixed-length DNA subsequence (here 18-22 bp, used as primer candidates) |
| Elite k-mer | A K-mer present in ≥ `design_min_conservation` fraction of target strains |
| Clean k-mer | Elite k-mer with zero hits in any background genome (strict mode) |
| Rescue mode | Adaptive fallback when strict filtering leaves too few clean k-mers |
| Anchor | A K-mer used as starting point for primer pair formation |
| Pigeonhole search | Fast exact-match acceleration for approximate string matching |
| IUPAC | Nucleotide ambiguity codes (R=AG, Y=CT, N=ACGT, etc.) |
| TaqMan probe | Fluorescent hydrolysis probe for qPCR (binds between primers) |
| Exclusivity panel | Background species used to ensure primers do not cross-react |
| Tm | Melting temperature of primer-DNA duplex |
| Amplicon | PCR product sequence |
| Cross-reactivity | Unwanted amplification of non-target DNA |
| In-silico PCR | Computational simulation of PCR amplification |
