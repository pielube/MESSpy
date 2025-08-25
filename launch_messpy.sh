#!/bin/bash

# MESSpy - Hybrid Plant Simulation Platform
# Launch script for Unix/Linux/macOS systems

echo ""
echo "========================================"
echo "   MESSpy Simulation Platform"
echo "   Hybrid Plant Analysis System"
echo "========================================"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Make the script executable
chmod +x wrapper.py

# Launch the application
python3 wrapper.py
