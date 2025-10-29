#!/bin/bash

# NVIDIA Protein Binding Design Workflow Launcher
# This script starts the multi-step workflow application

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  🧬 NVIDIA Protein Binding Design Workflow"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo "Please run setup.sh first:"
    echo "  bash setup.sh"
    exit 1
fi

# Activate virtual environment
source .venv/bin/activate

# Check if dependencies are installed
if ! python -c "import streamlit" 2>/dev/null; then
    echo "❌ Dependencies not installed!"
    echo "Installing dependencies..."
    pip install -r requirements.txt
fi

echo ""
echo ""
echo "🚀 Launching NVIDIA Binding Workflow Application..."
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  🌐 Application URL: http://localhost:8501"
echo "  🔬 Powered by NVIDIA AI Models"
echo ""
echo "  📊 Workflow Pipeline:"
echo "     ESMFold → RFDiffusion → ProteinMPNN → DiffDock"
echo ""
echo "  ⌨️  Press Ctrl+C to stop the application"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Set PYTHONPATH to project root
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Run Streamlit with the workflow app
streamlit run frontend/binding_workflow_app.py \
    --server.port=8501 \
    --server.address=localhost \
    --browser.gatherUsageStats=false \
    --theme.primaryColor="#1f77b4" \
    --theme.backgroundColor="#ffffff" \
    --theme.secondaryBackgroundColor="#f0f2f6"
