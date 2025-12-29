# Rational Primer Design Pipeline

### 📊 Project Growth
| Metric | Count |
| :--- | :--- |
| **Total Views** | 0|
| **Total Clones/Downloads** | 0|
## 🧬 About
**Rational Primer Design** is a comprehensive pipeline for designing PCR primers. It automates the selection process to ensure high specificity and efficiency for your biological workflows.

## 🚀 Installation

### Option 1: Quick Install (Windows)
1. Download the latest release.
2. Double-click **`INSTALL.bat`**.

### Option 2: Quick Install (Linux/Mac)
1. Open your terminal.
2. Run: `bash INSTALL.sh`

### Option 3: Developer Install (Python)
If you want to modify the code or install manually:

```bash
# Clone the repository
git clone [https://github.com/thanh727/rational-primer-design.git](https://github.com/thanh727/rational-primer-design.git)
cd rational-primer-design

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install .

💻 Usage
After installation, you can run the pipeline using the command line:

Bash

# Example command
python rational_design/main.py --input data/genes.fasta --output results/

(Note: Please refer to USER_MANUAL.txt for detailed parameters and options.)

📂 Project Structure
rational_design/: Main source code.

config/: Configuration files.

database/: Reference data.

traffic/: Traffic data logs (generated automatically).
