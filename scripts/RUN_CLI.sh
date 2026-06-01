#!/bin/bash

# 1. Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
cd "$ROOT_DIR"

# 2. Activate Virtual Environment automatically
if [ -d "venv" ]; then
    source venv/bin/activate
else
    echo "❌ Error: Virtual environment 'venv' not found."
    echo "   Please run ./INSTALL_UX.sh first."
    exit 1
fi

# 3. Pass all arguments ($@) to the Python CLI
#    This allows the user to type "./RUN_CLI.sh pipeline --out ..."
python3 -m rational_design.cli "$@"
