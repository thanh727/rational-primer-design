# Rational Primer Design Pipeline 🧬
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18103221.svg)](https://doi.org/10.5281/zenodo.18103221)

**A high-performance pipeline for the rational design and in‑silico validation of TaqMan qPCR assays.**

Rational Primer Design is a scalable, parallelized bioinformatics framework for **designing, optimizing, and validating TaqMan qPCR assays** across large-scale genomic datasets. The pipeline is engineered for reproducibility, computational efficiency, and applicability to pathogen surveillance, diagnostics, and translational research.

---

## 🚀 Key Features

- **Automated Primer & Probe Design**
  Identification of conserved genomic regions using a fast **2‑bit integer encoding** strategy for scalable and memory‑efficient sequence comparison.

- **In‑silico PCR Validation**
  High‑throughput screening of candidate assays against hundreds to thousands of background genomes using a *Turbo Pigeonhole*–based mismatch pruning algorithm.

- **High‑Performance Parallelization**
  Efficient utilization of all available CPU cores to accelerate candidate generation, filtering, and validation.

- **Optional Constraint Relaxation**
  Constraint relaxation is available only when explicitly enabled, so user-provided sensitivity and conservation thresholds are not silently weakened.

- **Deterministic & Reproducible**
  Deterministic sampling and execution guarantee identical results across repeated runs given identical inputs.

- **Cross‑Platform Support**
  Native execution on **Windows**, **macOS**, and **Linux** systems.

---

## 📦 Installation

### 🚀 Quick Start (Cloud Version)

Run the complete application directly in your browser using the provided **Google Colab notebook**. This option requires **no local installation** and is recommended for evaluation or rapid prototyping.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thanh727/rational-primer-design/blob/main/notebooks/Application_primer_design.ipynb)

---

### 🖥️ Desktop Installation (Local)

Download the latest release from the **Releases** tab of this repository, extract the archive, and follow the operating system–specific instructions below.

#### 🪟 Windows

1. Download and unzip the program archive.
2. Double‑click `INSTALL.bat` to:
   - Install all required dependencies
   - Register the `rational-design` command in the system environment

---

#### 🍎 macOS / 🐧 Linux

##### 1. Prerequisite: Tkinter

Tkinter is required for the graphical user interface.

**Ubuntu / Debian**
```bash
sudo apt-get update
sudo apt-get install python3-tk
```

**macOS**
If Python was installed via Homebrew, Tkinter is typically included. If missing:
```bash
brew install python-tk
```

##### 2. Install the pipeline
```bash
cd /path/to/rational_primer_design
chmod +x INSTALL_UX.sh RUN_APP_UX.sh
./INSTALL_UX.sh
```

---

## 🏃‍♂️ Running the Pipeline

### 🌐 Web App (Next.js + FastAPI)

The modern web interface runs the Python pipeline through a FastAPI backend and a Next.js frontend.

#### 🪟 Windows
Double-click:
```text
run_local.bat
```

#### 🍎 macOS / 🐧 Linux
```bash
chmod +x run_local.sh
./run_local.sh
```

Once started, open:
```text
http://127.0.0.1:3000
```

*Note: The web server automatically monitors active browser tabs. When you close all browser tabs, the server and all associated background processes will automatically shut down to save system resources.*

---

### 🖥️ Desktop (Interactive Terminal Wizard)

The interactive terminal wizard guides you step-by-step through configuring and running pipeline tasks.

#### 🪟 Windows
Double‑click:
```text
RUN_APP.bat
```

#### 🍎 macOS / 🐧 Linux
```bash
chmod +x RUN_APP_UX.sh
./RUN_APP_UX.sh
```

Alternatively, if installed as a package, run:
```bash
rational-design term
```

---

### ⌨️ Command‑Line Interface (CLI Mode)

If you installed the package, you can run the pipeline directly from your terminal using the globally registered `rational-design` command.

#### 1. Auto‑Download Mode (NCBI Fetch)

**Recommended for:** Starting from scratch and allowing the pipeline to automatically download genomes from NCBI.

Prepare two configuration files:
- `targets.json`
- `background.json`

**Example `targets.json`:**
```json
{
  "Salmonella": ["Salmonella enterica[Org] AND complete genome", 2.0, 50]
}
```

**Format:**
```
"Name": ["NCBI Query", MinGenomeSize_MB, MaxGenomeCount]
```

**Run the pipeline:**
```bash
rational-design pipeline --out "results_auto_test" --email "your_email@example.com" --target_config "targets.json" --bg_config "background.json"
```

> An email address is required by NCBI for genome downloads.

---

#### 2. Local Mode (Pre‑Downloaded Genomes)

**Recommended for:** Running on existing genome collections (e.g. servers or HPC clusters).

```bash
rational-design pipeline --out "results_local_test" --local_target "path/to/target_genomes_folder" --local_bg "path/to/background_genomes_folder"
```

**Optional:**
If a default configuration file (`config/parameters.json`) exists, no additional arguments are required.
To use custom settings:
```bash
--params "my_custom_settings.json"
```

---

## ⚙️ Configuration

All biological and computational parameters are configured interactively via the left control panel or JSON configuration files.

| Parameter | Default | Description |
|---------|---------|-------------|
| design_target_sampling_size | 0 | Number of target genomes used for design (0 = all genomes; maximum accuracy). |
| design_max_candidates | 10 | Number of primer–probe sets evaluated per design cycle. |
| min_sensitivity | 90.0 | Minimum percentage of target genomes detected. |
| primer_length_min / primer_length_max | 18 / 22 | Primer length search range (bp). |
| product_size_min | 100 | Minimum amplicon size (bp). |
| product_size_max | 400 | Maximum amplicon size (bp). |
| enable_blast | false | Annotate amplicons using NCBI BLAST when explicitly enabled (requires internet connection). |
| auto_relax_constraints | false | If true, retries with lower conservation/sensitivity thresholds after strict attempts fail. |

---

## 🛠 System Requirements

- **Python ≥ 3.9**
  Must be installed and accessible via the system `PATH`.

- **Node.js ≥ 20**
  Required for the Next.js frontend.

- **Internet connection (optional)**
  Required only for genome downloading and BLAST‑based annotation.

- **Hardware considerations**
  On a standard desktop or laptop equipped with:
  - **16 GB RAM**
  - **Intel Core i5‑9400 processor or equivalent**

  the pipeline can typically process:
  - Up to **150 target** and **150 background genomes** for **Gram‑positive bacteria**
  - Up to **50 target** and **50 background genomes** for **Gram‑negative bacteria**

  These limits scale proportionally with increased memory and CPU capacity.

---

## 📚 References

- Dieffenbach CW, Dveksler GS. *PCR Primer: A Laboratory Manual*. Cold Spring Harbor Laboratory Press.
- Thornton B, Basu C. Rapid and specific detection of bacteria by TaqMan PCR. *Journal of Microbiological Methods*.
- Untergasser A, et al. Primer3—new capabilities and interfaces. *Nucleic Acids Research*.

---

## 📚 Citation

If you use **Rational Primer Design** in your research, please cite the specific version used:

> Nguyen, T. (2025). *Rational Primer Design: High‑performance automated assay design pipeline* (Version 1.0.3). Zenodo. https://doi.org/10.5281/zenodo.18103221

**BibTeX:**
```bibtex
@software{Nguyen_Rational_Primer_Design_2025,
  author    = {Nguyen, Thanh},
  title     = {{Rational Primer Design}},
  version   = {1.0.3},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18103221},
  url       = {https://github.com/thanh727/rational-primer-design}
}
```
