# Rational Primer Design Pipeline 🧬

A high-performance, parallelized bioinformatics tool for designing and validating TaqMan PCR assays.

## 🚀 Features
* **Automated Design:** Mines conserved regions from target genomes using a 2-Bit integer algorithm.
* **In-Silico Validation:** Simulates PCR against hundreds of background genomes using "Turbo Pigeonhole" logic.
* **Parallel Processing:** Uses all available CPU cores for maximum speed.
* **Smart Optimization:** Automatically relaxes constraints if strict parameters fail to find candidates.
* **Reproducible:** Uses deterministic sampling so results are identical between runs.
* **Cross-Platform:** Runs natively on **Windows**, **macOS**, and **Linux**.

---

## 📦 Installation

### 🪟 Windows
1. Double-click **`INSTALL.bat`**.
   *(This installs Python dependencies and registers the `rational-design` command).*

### 🍎 macOS / 🐧 Linux
1. Open your **Terminal**.
2. Navigate to this folder:
   ```bash
   cd /path/to/rational_primer_design
   ```
3. Make the scripts executable (only needed once):
   ```bash
   chmod +x INSTALL.sh RUN_PIPELINE.sh
   ```
4. Run the installer:
   ```bash
   ./INSTALL.sh
   ```

## 🏃‍♂️ How to Run

### 🪟 Windows
Double-click **RUN_PIPELINE.bat**.

Follow the on-screen prompts to select your Project Name and Data Mode.

### 🍎 macOS / 🐧 Linux
Open Terminal and run:
```bash
./RUN_PIPELINE.sh
```
Follow the on-screen prompts.

## ⚙️ Configuration
You can customize the biological parameters by editing: `config/parameters.json`

| Parameter | Default | Description |
|---------|---------|-------------|
| design_target_sampling_size | 0 | Number of genomes to use for design. 0 = Use ALL (Highest Accuracy). |
| design_max_candidates | 10 | Number of primer pairs to attempt in the first cycle. |
| min_sensitivity | 95.0 | Minimum % of target genomes the primer must detect. |
| primer_length | 20 | Length of the primers (bp). |
| product_size_min | 100 | Minimum amplicon length. |
| enable_blast | true | Automatically annotate gene names via NCBI BLAST. |

## 📂 Folder Structure
- `rational_design/`: Source code (Python).
- `config/`: Configuration files.
- `database/`: Place your FASTA files here (target and background subfolders).
- `MyProject/`: Output folder (Generated after running).

## 🛠 Requirements
- Python 3.9+ (Must be installed and added to PATH).
- Internet Connection (Optional, only for downloading genomes from NCBI/BLAST).
