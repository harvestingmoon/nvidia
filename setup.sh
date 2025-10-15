#!/bin/bash

# Setup script for Protein Structure Prediction Application
# Run this script to set up the environment and install dependencies

echo "🧬 Setting up Protein Structure Prediction Application"
echo "=================================================="

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.8 or higher."
    exit 1
fi

echo "✅ Python 3 found: $(python3 --version)"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "❌ pip3 is not installed. Please install pip."
    exit 1
fi

echo "✅ pip3 found"

# Install requirements
echo "📦 Installing Python dependencies..."
pip3 install -r requirements.txt

if [ $? -eq 0 ]; then
    echo "✅ Dependencies installed successfully"
else
    echo "❌ Failed to install dependencies"
    exit 1
fi

# Create .env file from template if it doesn't exist
if [ ! -f ".env" ]; then
    echo "📝 Creating .env file from template..."
    cp .env.template .env
    echo ".env file created. Please edit it with your API credentials."
else
    echo " .env file already exists"
fi

echo ""
echo "🎉 Setup complete!"
echo ""
echo "Available Applications:"
echo "1. Enhanced App (Recommended): streamlit run app_v2.py --server.port 8502"
echo "2. Original App: streamlit run app.py"
echo "3. PDB Viewer: streamlit run pdb_viewer.py --server.port 8503"
echo ""
echo "Quick start options:"
echo "• Run: ./launch.sh (easiest - launches enhanced app)"
echo "• Or: streamlit run app_v2.py --server.port 8502"
echo ""
echo "URLs:"
echo "• Enhanced App: http://localhost:8502 (recommended)"
echo "• Original App: http://localhost:8501"
echo "• PDB Viewer: http://localhost:8503"
echo ""
echo "Happy protein folding! 🧬"
