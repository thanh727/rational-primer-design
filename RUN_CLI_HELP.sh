#!/bin/bash
# Activate the virtual environment if it exists
if [ -d "venv" ]; then
    source venv/bin/activate
fi

# Show the help menu for the CLI
python3 -m rational_design.cli --help