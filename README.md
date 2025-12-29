# Rational Primer Design Pipeline 🧬

A high-performance, parallelized bioinformatics tool for designing and validating TaqMan PCR assays.

### 📊 Project Growth
![Traffic Stats](traffic/traffic_badge.svg)

---

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
