#!/bin/bash

# NVIDIA Protein Structure Prediction App Launcher
# This script launches the Streamlit application for protein structure prediction

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  🧬 NVIDIA Protein Structure Prediction"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Activate virtual environment if it exists
if [ -d ".venv" ]; then
    echo "📦 Using virtual environment..."
    source .venv/bin/activate
fi

# Check if required packages are installed
echo "🔍 Checking dependencies..."
python -c "import streamlit, requests, py3Dmol" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "❌ Missing dependencies. Installing..."
    pip install -r requirements.txt
fi

# Launch the application
echo "🚀 Starting NVIDIA AI Application..."
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  🌐 Application URL: http://localhost:8502"
echo "  🔬 Powered by NVIDIA Cloud Functions"
echo ""
echo "  ⌨️  Press Ctrl+C to stop the application"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Set PYTHONPATH to project root
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

streamlit run frontend/app_v2.py --server.port 8502 --server.address localhost
