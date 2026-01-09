#!/bin/bash
# Launch the Binder Examples Gallery

cd "$(dirname "$0")/.."

echo "🧬 Starting NVIDIA Protein Binder Examples Gallery..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

streamlit run frontend/examples_gallery.py \
    --server.port=8502 \
    --server.address=localhost \
    --browser.gatherUsageStats=false
