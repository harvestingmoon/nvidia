#!/bin/bash

# Protein Structure Prediction App Launcher
# This script launches the Streamlit application for protein structure prediction

echo "🧬 Launching Protein Structure Prediction App"
echo "============================================"

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
echo "🚀 Starting Streamlit app..."
echo ""
echo "The app will be available at:"
echo "   → http://localhost:8502"
echo ""
echo "Press Ctrl+C to stop the application"
echo ""

streamlit run app_v2.py --server.port 8502 --server.address localhost
