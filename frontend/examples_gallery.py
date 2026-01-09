"""
Binder Examples Gallery
Showcase of successfully designed protein binders with 3D visualizations
"""

import streamlit as st
import sys
from pathlib import Path
import streamlit.components.v1 as components

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from frontend.app_v2 import create_3d_visualization

# Configure page
st.set_page_config(
    page_title="Binder Examples Gallery - NVIDIA",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# NVIDIA Theme CSS
st.markdown("""
<style>
    /* NVIDIA Brand Colors */
    :root {
        --nvidia-green: #76B900;
        --nvidia-dark: #1A1A1A;
    }
    
    h1, h2, h3 {
        color: #1A1A1A !important;
        font-family: 'NVIDIA Sans', Arial, sans-serif;
    }
    
    h1 {
        border-bottom: 3px solid #76B900;
        padding-bottom: 10px;
    }
    
    .stButton > button {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);
        color: white;
        border: none;
        border-radius: 6px;
        padding: 10px 24px;
        font-weight: 600;
    }
    
    .example-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
        border-left: 4px solid #76B900;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    
    .example-title {
        color: #1A1A1A;
        font-size: 24px;
        font-weight: 700;
        margin-bottom: 10px;
    }
    
    .example-description {
        color: #495057;
        font-size: 14px;
        line-height: 1.6;
        margin-bottom: 15px;
    }
    
    .metric-box {
        background: white;
        border-radius: 6px;
        padding: 10px;
        text-align: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    
    .metric-value {
        color: #76B900;
        font-size: 24px;
        font-weight: 700;
    }
    
    .metric-label {
        color: #6c757d;
        font-size: 12px;
        text-transform: uppercase;
    }
</style>
""", unsafe_allow_html=True)

# Example metadata
EXAMPLES = {
    "COVID_SPIKE_multimer_3.bin": {
        "name": "COVID-19 Spike Protein Binder",
        "target": "SARS-CoV-2 Spike Protein (RBD)",
        "description": "A computationally designed binder targeting the receptor-binding domain (RBD) of the SARS-CoV-2 spike protein. This binder could potentially neutralize the virus by blocking ACE2 receptor binding.",
        "use_case": "Antiviral Therapeutics",
        "confidence": "High",
        "applications": ["COVID-19 treatment", "Diagnostic tools", "Vaccine research"],
        "icon": "🦠"
    },
    "EGFR_multimer_1.bin": {
        "name": "EGFR Inhibitor Binder",
        "target": "Epidermal Growth Factor Receptor",
        "description": "A designed protein binder targeting EGFR, a receptor tyrosine kinase that plays a crucial role in cell proliferation. Overactive EGFR signaling is implicated in multiple cancers.",
        "use_case": "Cancer Therapeutics",
        "confidence": "Very High",
        "applications": ["Cancer treatment", "Targeted therapy", "Biomarker detection"],
        "icon": "🎗️"
    },
    "KRAS_G12D_multimer_1.bin": {
        "name": "KRAS G12D Mutant Binder",
        "target": "KRAS G12D Oncogenic Mutant",
        "description": "A precision binder designed to target the KRAS G12D mutation, one of the most common oncogenic drivers in pancreatic, lung, and colorectal cancers. This represents a 'undruggable' target made accessible.",
        "use_case": "Precision Oncology",
        "confidence": "Very High",
        "applications": ["Pancreatic cancer", "Lung cancer", "Colorectal cancer"],
        "icon": "🧬"
    },
    "pdl_1_multimer_3.bin": {
        "name": "PD-L1 Checkpoint Inhibitor",
        "target": "Programmed Death-Ligand 1 (PD-L1)",
        "description": "An immune checkpoint inhibitor targeting PD-L1, which cancer cells use to evade immune detection. Blocking PD-L1 can restore T-cell activity against tumors.",
        "use_case": "Cancer Immunotherapy",
        "confidence": "High",
        "applications": ["Immunotherapy", "Melanoma", "Lung cancer", "Bladder cancer"],
        "icon": "🛡️"
    },
    "5tpn_multimer_3.bin": {
        "name": "Insulin Receptor Binder",
        "target": "Insulin Receptor (5TPN)",
        "description": "A designed binder targeting the insulin receptor, with potential applications in diabetes research, biosensing, and development of novel insulin delivery systems.",
        "use_case": "Metabolic Disease Research",
        "confidence": "High",
        "applications": ["Diabetes research", "Biosensors", "Drug delivery"],
        "icon": "💉"
    }
}

def load_binary_pdb(file_path: Path) -> str:
    """Load and decode binary PDB file"""
    try:
        with open(file_path, 'rb') as f:
            pdb_bytes = f.read()
            pdb_content = pdb_bytes.decode('utf-8', errors='ignore')
        return pdb_content
    except Exception as e:
        st.error(f"Error loading {file_path.name}: {str(e)}")
        return None

def render_example_card(filename: str, metadata: dict, pdb_content: str):
    """Render a single example card with visualization"""
    st.markdown(f"""
    <div class="example-card">
        <div class="example-title">{metadata['icon']} {metadata['name']}</div>
        <div class="example-description">
            <strong>Target:</strong> {metadata['target']}<br>
            {metadata['description']}
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Metrics row
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown(f"""
        <div class="metric-box">
            <div class="metric-value">{metadata['use_case'].split()[0]}</div>
            <div class="metric-label">Use Case</div>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown(f"""
        <div class="metric-box">
            <div class="metric-value">{metadata['confidence']}</div>
            <div class="metric-label">Confidence</div>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown(f"""
        <div class="metric-box">
            <div class="metric-value">{len(metadata['applications'])}</div>
            <div class="metric-label">Applications</div>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        # Extract basic stats from PDB
        atom_count = pdb_content.count('ATOM')
        st.markdown(f"""
        <div class="metric-box">
            <div class="metric-value">{atom_count:,}</div>
            <div class="metric-label">Atoms</div>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # 3D Visualization
    with st.expander("🔬 View 3D Structure", expanded=True):
        if pdb_content:
            try:
                html_content = create_3d_visualization(pdb_content, color_by_plddt=True)
                components.html(html_content, height=600, scrolling=False)
            except Exception as e:
                st.error(f"Visualization error: {str(e)}")
        else:
            st.warning("Could not load structure data")
    
    # Applications
    with st.expander("📋 Potential Applications"):
        for app in metadata['applications']:
            st.markdown(f"• **{app}**")
    
    st.markdown("<hr style='margin: 30px 0; border: 1px solid #e9ecef;'>", unsafe_allow_html=True)

def main():
    # Header
    st.markdown("""
    <div style="text-align: center; padding: 20px 0;">
        <h1>🧬 PROTEIN BINDER EXAMPLES GALLERY</h1>
        <p style="font-size: 18px; color: #6c757d; margin-top: 10px;">
            Explore computationally designed protein binders targeting key therapeutic targets
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.title("🎯 Examples Gallery")
        st.markdown("---")
        
        st.markdown("""
        ### About These Examples
        
        These are **protein-protein complex structures** showing:
        
        - **Target protein** (disease-related)
        - **Designed binder** (therapeutic)
        - **Interface interactions**
        
        Each structure was predicted using **AlphaFold-Multimer** with confidence scores indicated by color.
        """)
        
        st.markdown("---")
        
        # Filter options
        st.subheader("Filter by Use Case")
        use_cases = list(set([ex['use_case'] for ex in EXAMPLES.values()]))
        selected_use_case = st.selectbox("Select category", ["All"] + use_cases)
        
        st.markdown("---")
        
        # Color legend
        st.subheader("pLDDT Color Legend")
        st.markdown("""
        <div style="background: #2d2d2d; padding: 15px; border-radius: 8px; color: white;">
            <div style="display: flex; align-items: center; margin-bottom: 8px;">
                <div style="width: 20px; height: 20px; background: #0053D6; margin-right: 10px; border-radius: 3px;"></div>
                <span>&gt;90: Very high confidence</span>
            </div>
            <div style="display: flex; align-items: center; margin-bottom: 8px;">
                <div style="width: 20px; height: 20px; background: #65CBF3; margin-right: 10px; border-radius: 3px;"></div>
                <span>70-90: Confident</span>
            </div>
            <div style="display: flex; align-items: center; margin-bottom: 8px;">
                <div style="width: 20px; height: 20px; background: #FFDB13; margin-right: 10px; border-radius: 3px;"></div>
                <span>50-70: Low confidence</span>
            </div>
            <div style="display: flex; align-items: center;">
                <div style="width: 20px; height: 20px; background: #FF7D45; margin-right: 10px; border-radius: 3px;"></div>
                <span>&lt;50: Very low</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Info banner
    st.info("💡 **Tip:** Click on any structure to interact - rotate, zoom, and inspect residues!")
    
    # Load examples directory
    examples_dir = Path(__file__).parent.parent / "examples_examples"
    
    if not examples_dir.exists():
        st.error(f"Examples directory not found: {examples_dir}")
        return
    
    # Render examples
    for filename, metadata in EXAMPLES.items():
        # Filter by use case
        if selected_use_case != "All" and metadata['use_case'] != selected_use_case:
            continue
        
        file_path = examples_dir / filename
        if file_path.exists():
            pdb_content = load_binary_pdb(file_path)
            if pdb_content:
                render_example_card(filename, metadata, pdb_content)
        else:
            st.warning(f"File not found: {filename}")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; color: #6c757d; padding: 20px;">
        <p>🚀 Powered by <strong style="color: #76B900;">NVIDIA NIM APIs</strong> | AlphaFold-Multimer | RFDiffusion | ProteinMPNN</p>
        <p style="font-size: 12px;">These are computational predictions for research purposes</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
