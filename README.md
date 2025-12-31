
# Rational Primer Design Pipeline 🧬
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18103221.svg)](https://doi.org/10.5281/zenodo.18103221)

**A high-performance pipeline for designing and validating TaqMan assays.**
...

A high-performance, parallelized bioinformatics framework for **designing, optimizing, and validating TaqMan qPCR assays** across large-scale genomic datasets.

---

## 🚀 Key Features

- **Automated Primer & Probe Design**  
  Identification of conserved genomic regions using a fast **2-bit integer encoding** strategy for scalable sequence comparison.

- **In-silico PCR Validation**  
  Screening of candidate assays against hundreds to thousands of background genomes using a *Turbo Pigeonhole*–based mismatch pruning algorithm.

- **High-Performance Parallelization**  
  Efficient utilization of all available CPU cores to accelerate candidate generation and validation.

- **Adaptive Constraint Relaxation**  
  Automatic relaxation of biological constraints when no valid primer–probe sets are found under strict parameter settings.

- **Deterministic & Reproducible**  
  Deterministic sampling ensures identical results across repeated runs given identical inputs.

- **Cross-Platform Support**  
  Native execution on **Windows**, **macOS**, and **Linux**.

---

## 📦 Installation

### 🚀 Quick Start (Cloud Version)

Run the full application directly in your browser using the provided Google Colab notebook (no local installation required):

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thanh727/rational-primer-design/blob/main/Application_primer_design.ipynb)

---

### 🖥️ Desktop Installation

Download the latest release from the **Releases** tab and follow the platform-specific instructions below.

#### 🪟 Windows

1. Download and unzip the program archive.
2. Double-click `INSTALL.bat` to install all dependencies and register the `rational-design` command.

---

#### 🍎 macOS / 🐧 Linux

##### 1. Prerequisite: Tkinter

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

##### 2. Install the pipeline
```bash
cd /path/to/rational_primer_design
chmod +x INSTALL_UX.sh RUN_APP_UX.sh
./INSTALL_UX.sh
```

---

## 🏃‍♂️ Running the Pipeline

### 🪟 Windows
Double-click:
```
RUN_APP.bat
```
Then follow the on-screen prompts.

### 🍎 macOS / 🐧 Linux
```bash
./RUN_APP_UX.sh
```
Follow the interactive prompts.

---

## ⚙️ Configuration

All biological and computational parameters are configured interactively via the left control panel:

| Parameter | Default | Description |
|----------|---------|-------------|
| design_target_sampling_size | 0 | Number of target genomes used for design (0 = all genomes; maximum accuracy). |
| design_max_candidates | 10 | Number of primer–probe sets evaluated per design cycle. |
| min_sensitivity | 90.0 | Minimum percentage of target genomes detected. |
| primer_length | 20 | Primer length (bp). |
| product_size_min | 100 | Minimum amplicon size (bp). |
| product_size_max | 400 | Maximum amplicon size (bp). |
| enable_blast | true | Annotate amplicons using NCBI BLAST (requires internet connection). |

---

## 🛠 System Requirements

- **Python ≥ 3.9**  
  Must be installed and available in the system `PATH`.

- **Internet connection (optional)**  
  Required only for genome downloading and BLAST-based annotation.

- **Hardware considerations**  
  Performance and genome capacity depend on available computational resources. On a standard desktop or laptop equipped with:
  - **16 GB RAM**
  - **Intel Core i5-9400 processor or equivalent**

  the pipeline can typically handle:
  - Up to **500 target genomes** and **500 background genomes** for **Gram-positive bacteria**
  - Up to **150 target genomes** and **150 background genomes** for **Gram-negative bacteria**

  These limits can be increased proportionally on systems with higher memory and computational capacity.

---

## 📚 References

- Dieffenbach CW, Dveksler GS. *PCR Primer: A Laboratory Manual*. Cold Spring Harbor Laboratory Press.
- Thornton B, Basu C. Rapid and specific detection of bacteria by TaqMan PCR. *Journal of Microbiological Methods*.
- Untergasser A, et al. Primer3—new capabilities and interfaces. *Nucleic Acids Research*.

## 📚 Citation
If you use Rational Primer Design in your research, please cite the specific version used:

> **Nguyen, T.** (2025). *Rational Primer Design: High-performance automated assay design pipeline (Version 1.0.3)*. 
> Zenodo. https://doi.org/10.5281/zenodo.18103221

Or use the BibTeX:
```bibtex
@software{Nguyen_Rational_Primer_Design_2025,
  author = {Nguyen, Thanh},
  title = {{Rational Primer Design}},
  version = {1.0.3},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18103221},
  url = {[https://github.com/thanh727/rational-primer-design](https://github.com/thanh727/rational-primer-design)}
}
