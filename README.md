# Rational Primer Design Pipeline 🧬

A high‑performance, parallelized bioinformatics framework for **designing, optimizing, and validating TaqMan qPCR assays** across large genomic datasets.

---

## 📊 Project Activity
![Traffic Stats](traffic/traffic_badge.svg)

---
## 🚀 Key Features

- **Automated Primer & Probe Design**  
  Identifies conserved genomic regions using a fast **2‑bit integer encoding** strategy for scalable sequence comparison.

- **In‑Silico PCR Validation**  
  Screens candidate assays against hundreds to thousands of background genomes using a *Turbo Pigeonhole*–based mismatch pruning algorithm.

- **High‑Performance Parallelization**  
  Fully utilizes available CPU cores to accelerate candidate generation and validation.

- **Adaptive Constraint Relaxation**  
  Automatically relaxes biological constraints when no valid candidates are found under strict parameters.

- **Deterministic & Reproducible**  
  Deterministic sampling ensures identical results across repeated runs with the same inputs.

- **Cross‑Platform Support**  
  Native execution on **Windows**, **macOS**, and **Linux**.

---

## 📦 Installation

---

### 🚀 Quick Start (Cloud Version)
Don't want to install anything? Run the full application in your browser using our free Google Colab notebook.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thanh727/rational-primer-design/blob/main/Application_primer_design.ipynb)

---

### 🖥️ Desktop Installation
Prefer to run locally on Windows/Linux? Download the latest release from the "Releases" tab.
### 🪟 Windows

1. Open **Command Prompt**.
2. Navigate to the project directory:
   ```bash
   cd /path/to/rational_primer_design
   ```
3. Install in editable mode:
   ```bash
   pip install -e .
   ```
   This installs all dependencies and registers the `rational-design` command.

---

### 🍎 macOS / 🐧 Linux

#### 1. Prerequisite: Tkinter

**Ubuntu / Debian**
```bash
sudo apt-get update
sudo apt-get install python3-tk
```

**macOS**
If Python was installed via Homebrew, Tkinter is usually included. If missing:
```bash
brew install python-tk
```

#### 2. Install the pipeline

```bash
cd /path/to/rational_primer_design
chmod +x INSTALL.sh RUN_PIPELINE.sh RUN_APP_UX.sh
./INSTALL.sh
```

---

## 🏃‍♂️ Running the Pipeline

### 🪟 Windows
Double‑click:
```
RUN_APP.bat
```

Follow the on‑screen prompts to select:
- Project name  
- Target/background genome mode  

### 🍎 macOS / 🐧 Linux
```bash
./RUN_APP_UX.sh
```
Follow the interactive prompts.

---

## ⚙️ Configuration

All biological and computational parameters are defined in:

```
config/parameters.json
```

| Parameter | Default | Description |
|---------|---------|-------------|
| design_target_sampling_size | 0 | Number of target genomes used for design (0 = all genomes; highest accuracy). |
| design_max_candidates | 10 | Number of primer–probe sets evaluated per design cycle. |
| min_sensitivity | 90.0 | Minimum percentage of target genomes detected. |
| primer_length | 20 | Primer length (bp). |
| product_size_min | 100 | Minimum amplicon size (bp). |
| product_size_max | 400 | Maximum amplicon size (bp). |
| enable_blast | true | Annotate amplicons using NCBI BLAST (requires internet). |

---

## 📂 Directory Structure

```
rational_design/     # Core Python source code
config/              # Configuration files
database/
 ├─ target/           # Target genome FASTA files
 └─ background/       # Background genome FASTA files
MyProject/            # Output directory (auto‑generated)
```

---

## 🛠 System Requirements

- **Python ≥ 3.9** (installed and available in PATH)
- Internet connection *(optional; required only for genome download and BLAST annotation)*

---

## 📚 Recommended References

- Dieffenbach CW, Dveksler GS. *PCR Primer: A Laboratory Manual*. Cold Spring Harbor Laboratory Press.  
- Thornton B, Basu C. Rapid and specific detection of bacteria by TaqMan PCR. *J Microbiol Methods*.  
- Untergasser A et al. Primer3—new capabilities and interfaces. *Nucleic Acids Research*.
