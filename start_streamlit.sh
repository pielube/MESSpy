#!/bin/bash

# Hybrid Plant Simulation Dashboard - Streamlit Launcher
# This script starts the Streamlit application for the hybrid plant simulation

echo "⚡ Starting Hybrid Plant Simulation Dashboard..."
echo "================================================"

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "❌ Error: Python is not installed or not in PATH"
    exit 1
fi

# Check if Streamlit is installed
if ! python -c "import streamlit" &> /dev/null; then
    echo "❌ Error: Streamlit is not installed"
    echo "Please install it with: pip install -r streamlit_requirements.txt"
    exit 1
fi

# Check if required files exist
if [ ! -f "streamlit_app.py" ]; then
    echo "❌ Error: streamlit_app.py not found"
    exit 1
fi

if [ ! -f "run_test_4.py" ]; then
    echo "❌ Error: run_test_4.py not found"
    exit 1
fi

# Check if configuration directory exists
if [ ! -d "input_test_4" ]; then
    echo "❌ Error: input_test_4 directory not found"
    exit 1
fi

echo "✅ All requirements met"
echo ""
echo "🚀 Launching Streamlit application..."
echo "📱 The application will open in your default web browser"
echo "🌐 URL: http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop the application"
echo ""

# Start Streamlit
streamlit run streamlit_app.py --server.headless false --server.port 8501
