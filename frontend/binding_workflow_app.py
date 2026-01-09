"""
Binding Protein Design Workflow Application
Multi-step workflow for designing and analyzing protein binders
AlphaFold2 → RFDiffusion → ProteinMPNN → AlphaFold Multimer
"""

import streamlit as st
import streamlit.components.v1 as components
import time
import json
import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
import plotly.graph_objects as go
import plotly.express as px
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from workflow.workflow_state import (
    WorkflowSession, WorkflowStage, StageStatus, 
    WorkflowValidator, TargetProteinData, BinderProteinData
)
from workflow.generative_pipeline import GenerativePipeline
from workflow.binding_analysis import (
    parse_pdb_content, find_interface_residues, 
    assess_binding_quality, combine_pdbs, generate_contact_map_data
)
from frontend.app_v2 import (
    call_nvidia_protein_api, validate_protein_sequence,
    create_3d_visualization, extract_pdb_from_response
)
from core.protein_models import PROTEIN_MODELS

# Configure page
st.set_page_config(
    page_title="NVIDIA Protein Binding Design Workflow",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# NVIDIA Theme Colors
NVIDIA_GREEN = "#76B900"
NVIDIA_DARK = "#1A1A1A"
NVIDIA_GRAY = "#666666"
NVIDIA_LIGHT_GRAY = "#E5E5E5"
NVIDIA_WHITE = "#FFFFFF"
NVIDIA_ACCENT = "#00D4AA"  # Teal accent
NVIDIA_WARNING = "#FFA500"


def check_has_plddt_scores(pdb_content: str) -> bool:
    """Check if PDB has meaningful pLDDT scores in B-factor column"""
    if not pdb_content:
        return False
    
    b_factors = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                b_factor = float(line[60:66].strip())
                b_factors.append(b_factor)
            except (ValueError, IndexError):
                continue
    
    if not b_factors:
        return False
    
    # Check if B-factors look like pLDDT scores (0-100 range with variation)
    min_b = min(b_factors)
    max_b = max(b_factors)
    avg_b = sum(b_factors) / len(b_factors)
    
    # pLDDT scores are typically 0-100 with meaningful variation
    has_variation = (max_b - min_b) > 5
    in_plddt_range = 0 <= min_b <= 100 and 0 <= max_b <= 100
    reasonable_avg = 20 <= avg_b <= 100
    
    return has_variation and in_plddt_range and reasonable_avg
NVIDIA_ERROR = "#FF3838"

# Custom CSS with NVIDIA Design System
st.markdown("""
<style>
    /* NVIDIA Brand Colors */
    :root {
        --nvidia-green: #76B900;
        --nvidia-dark: #1A1A1A;
        --nvidia-gray: #666666;
        --nvidia-light-gray: #E5E5E5;
        --nvidia-accent: #00D4AA;
    }
    
    /* Main App Styling */
    .stApp {
        background-color: #FFFFFF;
    }
    
    /* Headers */
    h1, h2, h3 {
        color: #1A1A1A !important;
        font-family: 'NVIDIA Sans', 'Helvetica Neue', Arial, sans-serif;
        font-weight: 600;
    }
    
    h1 {
        border-bottom: 3px solid #76B900;
        padding-bottom: 10px;
    }
    
    /* Progress Bar */
    .stProgress > div > div > div > div {
        background: linear-gradient(90deg, #76B900 0%, #00D4AA 100%);
    }
    
    /* Stage Status */
    .stage-complete {
        color: #76B900;
        font-weight: 700;
        font-size: 16px;
        text-shadow: 0 0 10px rgba(118, 185, 0, 0.3);
        padding: 10px;
        text-align: center;
        display: block;
        min-height: 50px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .stage-active {
        color: #00D4AA;
        font-weight: 700;
        font-size: 16px;
        animation: pulse 2s infinite;
        padding: 10px;
        text-align: center;
        display: block;
        min-height: 50px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .stage-pending {
        color: #666666;
        font-weight: 500;
        font-size: 16px;
        padding: 10px;
        text-align: center;
        display: block;
        min-height: 50px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .stage-failed {
        color: #FF3838;
        font-weight: 700;
        font-size: 16px;
        padding: 10px;
        text-align: center;
        display: block;
        min-height: 50px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    @keyframes pulse {
        0%, 100% { opacity: 1; }
        50% { opacity: 0.6; }
    }
    
    /* Metric Cards */
    .metric-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e5e5e5 100%);
        padding: 20px;
        border-radius: 8px;
        border-left: 4px solid #76B900;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        transition: transform 0.2s;
    }
    
    .metric-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(118, 185, 0, 0.2);
    }
    
    /* Buttons - Base styles */
    .stButton > button {
        border-radius: 12px !important;
        padding: 16px 8px !important;
        font-weight: 600 !important;
        font-size: 13px !important;
        line-height: 1.3 !important;
        min-height: 85px !important;
        height: 85px !important;
        transition: all 0.3s ease !important;
        white-space: pre-line !important;
        display: flex !important;
        flex-direction: column !important;
        align-items: center !important;
        justify-content: center !important;
        gap: 4px !important;
    }
    
    .stButton > button:hover:not(:disabled) {
        transform: translateY(-2px) !important;
    }
    
    /* Secondary Buttons */
    .stButton > button[kind="secondary"] {
        background: linear-gradient(135deg, #666666 0%, #4D4D4D 100%);
    }
    
    /* Radio Buttons */
    .stRadio > div {
        background-color: #f8f9fa;
        padding: 10px;
        border-radius: 6px;
        border: 1px solid #E5E5E5;
    }
    
    /* Fix radio button label colors - AGGRESSIVE targeting for all radio buttons */
    .stRadio label,
    .stRadio label *,
    .stRadio div[role="radiogroup"] label,
    .stRadio div[role="radiogroup"] label *,
    .stRadio span,
    div[data-baseweb="radio"] label,
    div[data-baseweb="radio"] label *,
    div[data-baseweb="radio"] span {
        color: #1A1A1A !important;
        font-weight: 500 !important;
    }
    
    .main .stRadio label {
        color: #1A1A1A !important;
        font-weight: 500 !important;
    }
    
    .main .stRadio > label > div[data-testid="stMarkdownContainer"] > p {
        color: #1A1A1A !important;
    }
    
    .main .stRadio div[role="radiogroup"] label {
        color: #1A1A1A !important;
    }
    
    .main .stRadio div[role="radiogroup"] label span {
        color: #1A1A1A !important;
    }
    
    /* Fix all text inputs and labels in MAIN area only */
    .main label, .main .stMarkdown p, .main .stMarkdown span {
        color: #1A1A1A !important;
    }
    
    /* Text areas in main area */
    .main .stTextArea label, .main .stTextInput label {
        color: #1A1A1A !important;
        font-weight: 500 !important;
    }
    
    /* Force all text in main area to be dark */
    .main div[data-testid="stMarkdownContainer"] p,
    .main div[data-testid="stMarkdownContainer"] span {
        color: #1A1A1A !important;
    }
    
    /* Keep sidebar text WHITE/LIGHT */
    section[data-testid="stSidebar"] label,
    section[data-testid="stSidebar"] .stMarkdown p,
    section[data-testid="stSidebar"] .stMarkdown span,
    section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] p {
        color: #E5E5E5 !important;
    }
    
    /* Sidebar headers */
    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3 {
        color: #FFFFFF !important;
    }
    
    /* Sidebar warning text */
    section[data-testid="stSidebar"] .stAlert {
        color: #1A1A1A !important;
    }
    
    /* Keep sidebar RADIO buttons BLACK (visible on light background) */
    section[data-testid="stSidebar"] .stRadio div[role="radiogroup"] label,
    section[data-testid="stSidebar"] .stRadio div[role="radiogroup"] label span {
        color: #1A1A1A !important;
    }
    
    .stRadio > div [role="radiogroup"] label {
        color: #1A1A1A;
        font-weight: 500;
    }
    
    /* Text Inputs */
    .stTextInput > div > div > input,
    .stTextArea textarea {
        border: 2px solid #E5E5E5;
        border-radius: 6px;
        padding: 10px;
        font-size: 14px;
        transition: border-color 0.3s;
    }
    
    .stTextInput > div > div > input:focus,
    .stTextArea textarea:focus {
        border-color: #76B900;
        box-shadow: 0 0 0 2px rgba(118, 185, 0, 0.1);
    }
    
    /* Sidebar */
    .css-1d391kg, [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1A1A1A 0%, #2D2D2D 100%);
    }
    
    /* Sidebar buttons - consistent NVIDIA green styling */
    section[data-testid="stSidebar"] .stButton > button {
        background: linear-gradient(135deg, #76B900 0%, #5a8f00 100%) !important;
        color: white !important;
        border: none !important;
        min-height: 45px !important;
        height: auto !important;
        padding: 10px 16px !important;
    }
    
    section[data-testid="stSidebar"] .stButton > button:hover {
        background: linear-gradient(135deg, #8ed100 0%, #76B900 100%) !important;
    }
    
    .css-1d391kg h1, .css-1d391kg h2, .css-1d391kg h3,
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 {
        color: #76B900 !important;
    }
    
    .css-1d391kg p, .css-1d391kg label,
    [data-testid="stSidebar"] p, [data-testid="stSidebar"] label {
        color: #E5E5E5 !important;
    }
    
    /* Success Messages */
    .stSuccess {
        background-color: rgba(118, 185, 0, 0.1);
        border-left: 4px solid #76B900;
        color: #1A1A1A;
    }
    
    /* Info Messages */
    .stInfo {
        background-color: rgba(0, 212, 170, 0.1);
        border-left: 4px solid #00D4AA;
        color: #1A1A1A;
    }
    
    /* Warning Messages */
    .stWarning {
        background-color: rgba(255, 165, 0, 0.1);
        border-left: 4px solid #FFA500;
        color: #1A1A1A;
    }
    
    /* Error Messages */
    .stError {
        background-color: rgba(255, 56, 56, 0.1);
        border-left: 4px solid #FF3838;
        color: #1A1A1A;
    }
    
    /* Expander */
    .streamlit-expanderHeader {
        background-color: #f8f9fa;
        border-radius: 6px;
        border: 1px solid #E5E5E5;
        font-weight: 600;
        color: #1A1A1A;
    }
    
    .streamlit-expanderHeader:hover {
        border-color: #76B900;
    }
    
    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background-color: #f8f9fa;
        padding: 8px;
        border-radius: 8px;
    }
    
    .stTabs [data-baseweb="tab"] {
        background-color: transparent;
        border: 2px solid #E5E5E5;
        border-radius: 6px;
        color: #666666;
        font-weight: 600;
        padding: 8px 16px;
    }
    
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);
        border-color: #76B900;
        color: white;
    }
    
    /* Metrics */
    [data-testid="stMetricValue"] {
        color: #76B900;
        font-size: 32px;
        font-weight: 700;
    }
    
    [data-testid="stMetricLabel"] {
        color: #666666;
        font-weight: 600;
        font-size: 14px;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    /* Download Button */
    .stDownloadButton > button {
        background: linear-gradient(135deg, #00D4AA 0%, #00A88A 100%);
        color: white;
        border: none;
        font-weight: 600;
    }
    
    .stDownloadButton > button:hover {
        background: linear-gradient(135deg, #00A88A 0%, #00D4AA 100%);
    }
    
    /* Divider */
    hr {
        border: none;
        height: 2px;
        background: linear-gradient(90deg, #76B900 0%, #00D4AA 50%, #76B900 100%);
        margin: 30px 0;
    }
    
    /* File Uploader */
    [data-testid="stFileUploader"] {
        border: 2px dashed #76B900;
        border-radius: 8px;
        padding: 20px;
        background-color: rgba(118, 185, 0, 0.05);
    }
    
    /* Checkbox */
    .stCheckbox > label {
        color: #1A1A1A;
        font-weight: 500;
    }
    
    /* Select Box */
    .stSelectbox > div > div {
        border: 2px solid #E5E5E5;
        border-radius: 6px;
    }
    
    .stSelectbox > div > div:focus-within {
        border-color: #76B900;
        box-shadow: 0 0 0 2px rgba(118, 185, 0, 0.1);
    }
    
    /* Slider */
    .stSlider > div > div > div {
        background-color: #76B900;
    }
    
    /* NVIDIA Branding */
    .nvidia-header {
        background: linear-gradient(135deg, #1A1A1A 0%, #2D2D2D 100%);
        padding: 20px;
        border-radius: 8px;
        margin-bottom: 30px;
        border-left: 6px solid #76B900;
    }
    
    .nvidia-header h1 {
        color: #FFFFFF !important;
        border: none !important;
        margin: 0;
    }
    
    .nvidia-header p {
        color: #76B900;
        font-weight: 600;
        margin: 10px 0 0 0;
    }
    
    /* Quality Score Badge */
    .quality-badge {
        display: inline-block;
        padding: 8px 16px;
        border-radius: 20px;
        font-weight: 700;
        font-size: 18px;
        text-align: center;
    }
    
    .quality-excellent {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);
        color: white;
    }
    
    .quality-good {
        background: linear-gradient(135deg, #00D4AA 0%, #00A88A 100%);
        color: white;
    }
    
    .quality-moderate {
        background: linear-gradient(135deg, #FFA500 0%, #CC8400 100%);
        color: white;
    }
    
    .quality-poor {
        background: linear-gradient(135deg, #FF3838 0%, #CC2D2D 100%);
        color: white;
    }
</style>
""", unsafe_allow_html=True)


def render_pipeline_diagram():
    """Render the pipeline diagram using native Streamlit components"""
    
    # Custom CSS for pipeline styling
    st.markdown("""
    <style>
        .pipeline-header {
            text-align: center;
            padding: 20px 0;
        }
        .pipeline-header h1 {
            color: #76B900 !important;
            font-size: 2.5rem;
            font-weight: 900;
            letter-spacing: -1px;
            border: none !important;
        }
        .pipeline-subtitle {
            color: #888;
            font-size: 1.1rem;
        }
        .pipeline-divider {
            height: 2px;
            background: linear-gradient(to right, transparent, #76B900, transparent);
            margin: 15px auto;
            width: 60%;
        }
        .step-card {
            background: linear-gradient(135deg, #1a1a1a 0%, #111 100%);
            border: 2px solid #333;
            border-radius: 12px;
            padding: 0;
            transition: all 0.3s ease;
            height: 100%;
        }
        .step-card:hover {
            border-color: #76B900;
            box-shadow: 0 0 30px rgba(118, 185, 0, 0.3);
        }
        .step-header {
            background: linear-gradient(135deg, #222 0%, #1a1a1a 100%);
            padding: 15px;
            border-radius: 10px 10px 0 0;
            border-bottom: 2px solid #333;
        }
        .step-number {
            display: inline-block;
            background: rgba(118, 185, 0, 0.15);
            color: #76B900;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.75rem;
            font-weight: bold;
            margin-bottom: 8px;
        }
        .step-title {
            color: #fff;
            font-size: 1.1rem;
            font-weight: 700;
        }
        .step-content {
            padding: 15px;
            color: #ccc;
        }
        .model-tag {
            display: inline-block;
            background: rgba(118, 185, 0, 0.1);
            color: #76B900;
            padding: 5px 10px;
            border-radius: 4px;
            font-size: 0.8rem;
            margin: 3px;
            border: 1px solid rgba(118, 185, 0, 0.3);
        }
        .arrow-col {
            display: flex;
            align-items: center;
            justify-content: center;
            color: #76B900;
            font-size: 2rem;
        }
        .legend-box {
            display: inline-block;
            width: 18px;
            height: 18px;
            border-radius: 3px;
            margin-right: 8px;
            vertical-align: middle;
        }
        .app-card {
            background: linear-gradient(135deg, #1a1a1a 0%, #111 100%);
            border: 1px solid #333;
            border-radius: 12px;
            padding: 20px;
            height: 100%;
        }
        .app-card:hover {
            border-color: #76B900;
        }
        .app-title {
            color: #76B900;
            font-size: 1.1rem;
            font-weight: bold;
            margin-bottom: 10px;
        }
        .app-desc {
            color: #999;
            font-size: 0.9rem;
            line-height: 1.5;
        }
    </style>
    """, unsafe_allow_html=True)
    
    # Header
    st.markdown("""
    <div class="pipeline-header">
        <h1>🧬 PROTEIN BINDER DESIGN</h1>
        <div class="pipeline-divider"></div>
        <p class="pipeline-subtitle">AI-Powered Computational Pipeline</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Legend
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown("""
        <div style="text-align: center; color: #666; margin-bottom: 20px;">
            <span><span class="legend-box" style="background: #76B900;"></span><strong>Green:</strong> Target Protein</span>
            &nbsp;&nbsp;&nbsp;&nbsp;
            <span><span class="legend-box" style="background: #FF6600;"></span><strong>Orange:</strong> Designed Binder</span>
        </div>
        """, unsafe_allow_html=True)
    
    # Pipeline Steps Definition
    pipeline_steps = [
        {
            "num": 1,
            "title": "Target Input",
            "stage": WorkflowStage.TARGET_INPUT,
            "models": ["Amino Acid Sequence", "PDB File", "PDB ID"],
            "description": "Enter an amino acid sequence, upload a PDB file, or specify a PDB ID. The app validates your input and lets you specify binding site residues.",
            "pdb": None,  # No 3D view for input
            "color": "#76B900"
        },
        {
            "num": 2,
            "title": "Target Structure Prediction",
            "stage": WorkflowStage.TARGET_PREDICTION,
            "models": ["AlphaFold3", "AlphaFold2", "OpenFold3"],
            "description": "Predict the 3D structure of your target protein. Get visual representation and confidence scores.",
            "pdb": "1CRN",
            "color": "#76B900"
        },
        {
            "num": 3,
            "title": "Binder Scaffold Design",
            "stage": WorkflowStage.BINDER_SCAFFOLD_DESIGN,
            "models": ["RFDiffusion"],
            "description": "Design a protein backbone (scaffold) that binds your target. Predicts the 3D shape without sequence information.",
            "pdb": "1L2Y",
            "color": "#FF6600"
        },
        {
            "num": 4,
            "title": "Sequence Design",
            "stage": WorkflowStage.BINDER_SEQUENCE_DESIGN,
            "models": ["ProteinMPNN"],
            "description": "Takes the backbone structure and predicts possible amino acid sequences for the binder protein.",
            "pdb": "1L2Y",
            "color": "#FF6600"
        },
        {
            "num": 5,
            "title": "Complex Prediction",
            "stage": WorkflowStage.COMPLEX_PREDICTION,
            "models": ["AlphaFold2 Multimer"],
            "description": "Combines target and binder structures to analyze their interaction. Get visual representation and confidence scores.",
            "pdb": "1IGT",
            "color": "#76B900"
        }
    ]
    
    st.markdown("### 🔬 Pipeline Steps")
    st.markdown("**Click on any step to navigate directly to that stage:**")
    st.markdown("---")
    
    # Render steps in a row
    cols = st.columns([3, 1, 3, 1, 3, 1, 3, 1, 3])
    
    for i, step in enumerate(pipeline_steps):
        col_idx = i * 2  # 0, 2, 4, 6, 8
        
        with cols[col_idx]:
            # Create clickable card using a button
            st.markdown(f"""
            <div class="step-card">
                <div class="step-header">
                    <div class="step-number">STEP {step['num']}</div>
                    <div class="step-title">{step['title']}</div>
                </div>
                <div class="step-content">
                    <div style="margin-bottom: 10px;">
                        {''.join([f'<span class="model-tag">{m}</span>' for m in step['models']])}
                    </div>
                    <p style="font-size: 0.85rem; color: #aaa;">{step['description'][:100]}...</p>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Navigation button
            if st.button(f"Go to Step {step['num']} →", key=f"pipeline_step_{step['num']}", use_container_width=True):
                st.session_state.workflow_session.advance_to_stage(step['stage'])
                st.session_state.show_pipeline = False
                st.rerun()
        
        # Add arrow between steps (except after last step)
        if i < len(pipeline_steps) - 1:
            with cols[col_idx + 1]:
                st.markdown("""
                <div style="display: flex; align-items: center; justify-content: center; height: 200px; color: #76B900; font-size: 2rem;">
                    →
                </div>
                """, unsafe_allow_html=True)
    
    st.markdown("---")
    
    # Applications Section
    st.markdown("### 🎯 Why Protein Binding Matters")
    st.markdown("Understanding the real-world impact of protein binder design")
    
    app_cols = st.columns(3)
    
    applications = [
        {
            "title": "🚫 Blocking & Inhibition",
            "desc": "**Most Common Use:** Stop harmful proteins\n\n• Cancer: Block growth signals\n• COVID-19: Bind spike protein\n• Autoimmune: Block inflammation\n\n*Example: Humira ($20B+ sales)*"
        },
        {
            "title": "⚡ Activation & Enhancement", 
            "desc": "**Purpose:** Make proteins work better\n\n• Enzyme activation\n• Stabilize therapeutics\n• Trigger signaling pathways"
        },
        {
            "title": "🎯 Delivery & Targeting",
            "desc": "**Purpose:** Guide drugs to locations\n\n• Target cancer cells directly\n• Cross blood-brain barrier\n• Organ-specific delivery"
        },
        {
            "title": "🔬 Detection & Diagnostics",
            "desc": "**Purpose:** Find specific proteins\n\n• COVID rapid tests\n• Cancer screening\n• Pregnancy tests"
        },
        {
            "title": "🔧 Modification & Control",
            "desc": "**Purpose:** Change protein behavior\n\n• Extend drug half-life\n• Control localization\n• Regulate interactions"
        },
        {
            "title": "✨ Success Story: Insulin",
            "desc": "**Design binders that:**\n\n• Make insulin last longer\n• Prevent degradation\n• Control receptor binding\n\n*Design in minutes vs months!*"
        }
    ]
    
    for i, app in enumerate(applications):
        with app_cols[i % 3]:
            st.markdown(f"""
            <div class="app-card">
                <div class="app-title">{app['title']}</div>
                <div class="app-desc">{app['desc'].replace(chr(10), '<br>')}</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown("")  # Spacing
    
    st.markdown("---")
    
    # Back button
    col1, col2, col3 = st.columns([1, 1, 1])
    with col2:
        if st.button("🔙 Back to Workflow", use_container_width=True, type="primary"):
            st.session_state.show_pipeline = False
            st.rerun()


def initialize_session_state():
    """Initialize session state variables"""
    if 'workflow_session' not in st.session_state:
        st.session_state.workflow_session = WorkflowSession.create_new()
    
    if 'api_key' not in st.session_state:
        st.session_state.api_key = None
    
    if 'demo_mode' not in st.session_state:
        st.session_state.demo_mode = False
    
    # Show pipeline diagram by default on first load
    if 'show_pipeline' not in st.session_state:
        st.session_state.show_pipeline = True
    
    # Examples gallery state
    if 'show_examples' not in st.session_state:
        st.session_state.show_examples = False


def render_progress_stepper():
    """Render workflow progress stepper with clickable stages"""
    session = st.session_state.workflow_session
    
    stages = [
        ("Target Input", WorkflowStage.TARGET_INPUT),
        ("Target Structure", WorkflowStage.TARGET_PREDICTION),
        ("Binder Scaffold", WorkflowStage.BINDER_SCAFFOLD_DESIGN),
        ("Sequence Design", WorkflowStage.BINDER_SEQUENCE_DESIGN),
        ("Complex Prediction", WorkflowStage.COMPLEX_PREDICTION),
        ("Results", WorkflowStage.RESULTS)
    ]
    
    st.markdown("### Workflow Progress")
    
    # Custom CSS for stage-specific button styling
    stage_styles = """
    <style>
    /* Normalize button container spacing */
    div[data-testid="column"] > div > div[data-testid="stVerticalBlock"] > div {
        gap: 0 !important;
    }
    
    /* Remove form padding for clickable buttons */
    div[data-testid="column"] form {
        margin: 0 !important;
        padding: 0 !important;
    }
    
    /* Base styling for all stage buttons */
    div[data-testid="column"] button {
        border-radius: 12px !important;
        padding: 12px 6px !important;
        font-weight: 600 !important;
        font-size: 12px !important;
        line-height: 1.3 !important;
        min-height: 95px !important;
        height: 95px !important;
        max-height: 95px !important;
        width: 100% !important;
        margin: 0 !important;
        transition: all 0.3s ease !important;
        white-space: normal !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        display: flex !important;
        flex-direction: column !important;
        align-items: center !important;
        justify-content: center !important;
        word-wrap: break-word !important;
        hyphens: auto !important;
    }
    
    /* Completed stages - NVIDIA Green */
    button[kind="secondary"]:not(:disabled) {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%) !important;
        color: white !important;
        border: 2px solid #76B900 !important;
        box-shadow: 0 4px 8px rgba(118, 185, 0, 0.3) !important;
    }
    button[kind="secondary"]:not(:disabled):hover {
        background: linear-gradient(135deg, #5A8F00 0%, #76B900 100%) !important;
        box-shadow: 0 6px 12px rgba(118, 185, 0, 0.5) !important;
        transform: translateY(-2px) !important;
    }
    
    /* In-progress stages - NVIDIA Green (lighter shade) */
    button[kind="primary"]:disabled {
        background: linear-gradient(135deg, #8CC63F 0%, #76B900 100%) !important;
        color: white !important;
        border: 2px solid #8CC63F !important;
        box-shadow: 0 4px 8px rgba(140, 198, 63, 0.4) !important;
        opacity: 1 !important;
        animation: pulse-green 2s ease-in-out infinite !important;
    }
    
    @keyframes pulse-green {
        0%, 100% { 
            box-shadow: 0 4px 8px rgba(140, 198, 63, 0.4) !important;
        }
        50% { 
            box-shadow: 0 6px 16px rgba(118, 185, 0, 0.6) !important;
        }
    }
    
    /* Pending stages - Gray */
    button[kind="secondary"]:disabled {
        background: linear-gradient(135deg, #E5E5E5 0%, #CCCCCC 100%) !important;
        color: #666666 !important;
        border: 2px solid #E5E5E5 !important;
        opacity: 0.7 !important;
        box-shadow: none !important;
    }
    </style>
    """
    st.markdown(stage_styles, unsafe_allow_html=True)
    
    # Use columns with custom styled buttons
    cols = st.columns(len(stages))
    
    for idx, (col, (label, stage)) in enumerate(zip(cols, stages)):
        with col:
            status = session.stage_statuses.get(stage.value, "not_started")
            
            # Determine icon and button type
            # Priority: COMPLETED status takes precedence over current stage
            if status == StageStatus.COMPLETED.value:
                icon = "✅"
                clickable = True
                button_type = "secondary"
            elif status == StageStatus.IN_PROGRESS.value:
                icon = "🔄"
                clickable = False
                button_type = "primary"
            elif status == StageStatus.FAILED.value:
                icon = "❌"
                clickable = False
                button_type = "secondary"
            else:
                icon = "⭕"
                clickable = False
                button_type = "secondary"
            
            # Create button with icon / label format
            button_label = f"{icon} {label}"
            
            if clickable:
                # Clickable button for completed stages
                if st.button(
                    button_label,
                    key=f"nav_{stage.value}",
                    use_container_width=True,
                    type=button_type,
                    help="Click to view this stage"
                ):
                    session.advance_to_stage(stage)
                    st.rerun()
            else:
                # Non-clickable display for other stages
                st.button(
                    button_label,
                    key=f"display_{stage.value}",
                    use_container_width=True,
                    disabled=True,
                    type=button_type
                )
    
    st.progress(calculate_overall_progress(session))



def calculate_overall_progress(session: WorkflowSession) -> float:
    """Calculate overall workflow progress"""
    total_stages = len(WorkflowStage)
    completed = sum(1 for status in session.stage_statuses.values() 
                   if status == StageStatus.COMPLETED.value)
    return completed / total_stages


def render_sidebar():
    """Render sidebar with session management"""
    session = st.session_state.workflow_session
    
    # NVIDIA Logo in Sidebar - Clickable to return to start
    try:
        # Display logo first
        st.sidebar.image("image/nvidia.jpg", use_container_width=True)
    except Exception as e:
        st.sidebar.markdown("# NVIDIA")
    
    # Prediction Page button - only show when pipeline is visible
    if st.session_state.get('show_pipeline', False):
        if st.sidebar.button("Prediction Page", key="logo_home", use_container_width=True, help="Click to go to Prediction Page"):
            st.session_state.show_pipeline = False
            st.session_state._stage_nav_clicked = True  # Signal to close gallery
            session.advance_to_stage(WorkflowStage.TARGET_INPUT)
            st.rerun()
    else:
        # View Pipeline Diagram button - only on prediction page
        pipeline_path = Path(__file__).parent / "pipeline.html"
        if pipeline_path.exists():
            if st.sidebar.button("View Pipeline Diagram", key="pipeline_view", use_container_width=True, help="Click to view interactive pipeline diagram"):
                st.session_state.show_pipeline = True
                st.rerun()
    
    st.sidebar.title("Binding Workflow")
    
    # Project info
    st.sidebar.subheader("Project")
    project_name = st.sidebar.text_input(
        "Project Name",
        value=session.project_name,
        key="project_name_input"
    )
    if project_name != session.project_name:
        session.project_name = project_name
        session.last_updated = datetime.now().isoformat()
    
    st.sidebar.caption(f"Session ID: {session.session_id}")
    st.sidebar.caption(f"Created: {session.created_at[:19]}")
    st.sidebar.caption(f"Updated: {session.last_updated[:19]}")
    
    # NVIDIA-styled divider
    st.sidebar.markdown("""
        <hr style="height:3px;border:none;background:linear-gradient(90deg, #76B900 0%, #00D4AA 100%);margin:20px 0;" />
    """, unsafe_allow_html=True)
    
    # API Configuration
    st.sidebar.subheader("Configuration")
    
    demo_mode = st.sidebar.checkbox(
        "Demo Mode (No API)",
        value=st.session_state.demo_mode,
        help="Use demo mode to test the interface without API calls"
    )
    st.session_state.demo_mode = demo_mode
    
    if not demo_mode:
        # Get API key from environment variables (loaded from .env)
        default_api_key = os.getenv("NGC_CLI_API_KEY") or os.getenv("NVIDIA_API_KEY") or ""
        
        api_key = st.sidebar.text_input(
            "NVIDIA API Key",
            type="password",
            value=st.session_state.api_key or default_api_key,
            help="Get your API key from https://build.nvidia.com/"
        )
        if api_key:
            st.session_state.api_key = api_key
        
        # If no API key set yet, use default from env
        if not st.session_state.api_key:
            st.session_state.api_key = default_api_key
            
        # Show warning if no API key is configured
        if not st.session_state.api_key:
            st.sidebar.warning("⚠️ No API key configured. Please add NVIDIA_API_KEY to your .env file or enter it above.")
    
    # NVIDIA-styled divider
    st.sidebar.markdown("""
        <hr style="height:3px;border:none;background:linear-gradient(90deg, #76B900 0%, #00D4AA 100%);margin:20px 0;" />
    """, unsafe_allow_html=True)
    
    # Examples Gallery - toggle button
    st.sidebar.subheader("📚 Examples Gallery")
    if st.session_state.get('show_examples', False):
        # Show "Prediction Page" button when viewing examples
        if st.sidebar.button("🔙 Prediction Page", use_container_width=True, help="Return to prediction workflow"):
            st.session_state.show_examples = False
            st.session_state.show_pipeline = False  # Ensure pipeline view is closed
            st.session_state._stage_nav_clicked = True
            session.advance_to_stage(WorkflowStage.TARGET_INPUT)
            st.rerun()
    else:
        # Show "View Binder Examples" button when in prediction page
        if st.sidebar.button("🧬 View Binder Examples", use_container_width=True, help="View gallery of successfully designed binders"):
            st.session_state.show_examples = True
            st.rerun()


def save_session():
    """Save current session"""
    session = st.session_state.workflow_session
    filename = f"workflow_{session.session_id}.json"
    
    try:
        with open(filename, 'w') as f:
            f.write(session.to_json())
        st.sidebar.success(f"✅ Saved: {filename}")
    except Exception as e:
        st.sidebar.error(f"❌ Save failed: {str(e)}")


def load_session():
    """Load session from file"""
    uploaded_file = st.sidebar.file_uploader("Upload session JSON", type=['json'])
    
    if uploaded_file:
        try:
            json_content = uploaded_file.read().decode('utf-8')
            session = WorkflowSession.from_json(json_content)
            st.session_state.workflow_session = session
            st.sidebar.success("✅ Session loaded!")
            st.rerun()
        except Exception as e:
            st.sidebar.error(f"❌ Load failed: {str(e)}")


def load_pipeline_results(session: WorkflowSession, folder_path: str) -> bool:
    """
    Load pipeline results from an existing output folder.
    
    Args:
        session: WorkflowSession to populate
        folder_path: Path to output folder (e.g., 5tpn_AF2_output)
    
    Returns:
        True if loading succeeded, False otherwise
    """
    import os
    import glob
    
    try:
        # Check if folder exists
        if not os.path.exists(folder_path):
            st.error(f"Folder not found: {folder_path}")
            return False
        
        # Extract project name from folder
        folder_name = os.path.basename(folder_path.rstrip('/'))
        if '_AF2_output' in folder_name:
            project_name = folder_name.replace('_AF2_output', '')
            model_used = 'AlphaFold2'
        elif '_OF3_output' in folder_name:
            project_name = folder_name.replace('_OF3_output', '')
            model_used = 'OpenFold3'
        elif '_AF3_output' in folder_name:
            project_name = folder_name.replace('_AF3_output', '')
            model_used = 'AlphaFold3'
        else:
            # Try to infer from files in folder
            project_name = folder_name.replace('_output', '')
            # Check file names for model hints
            files_in_folder = os.listdir(folder_path)
            files_str = ' '.join(files_in_folder).lower()
            if 'af2' in files_str or 'alphafold2' in files_str:
                model_used = 'AlphaFold2'
            elif 'of3' in files_str or 'openfold' in files_str:
                model_used = 'OpenFold3'
            elif 'af3' in files_str or 'alphafold3' in files_str:
                model_used = 'AlphaFold3'
            else:
                model_used = 'Uploaded PDB'
        
        # Load target structure (Step 1)
        target_pdb = os.path.join(folder_path, f"{project_name}_first_structure.pdb")
        if os.path.exists(target_pdb):
            with open(target_pdb, 'r') as f:
                session.target.pdb_content = f.read()
            session.target.structure_predicted = True
            session.target.structure_file_path = target_pdb
            session.target.model_used = model_used
            
            # Extract sequence from PDB
            from workflow.binding_analysis import parse_pdb_content
            atoms = parse_pdb_content(session.target.pdb_content)
            if atoms:
                # Get unique residues
                residues = {}
                for atom in atoms:
                    chain = atom.chain
                    res_num = atom.residue_number
                    res_name = atom.residue_name
                    key = (chain, res_num)
                    if key not in residues:
                        residues[key] = res_name
                
                # Convert 3-letter codes to 1-letter
                aa_map = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }
                sorted_residues = sorted(residues.items(), key=lambda x: x[0][1])
                sequence = ''.join([aa_map.get(res[1], 'X') for res in sorted_residues])
                session.target.sequence = sequence
            
            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
        
        # Load binder scaffold (Step 2)
        scaffold_pdb = os.path.join(folder_path, f"{project_name}_RFD_prediction.pdb")
        if os.path.exists(scaffold_pdb):
            with open(scaffold_pdb, 'r') as f:
                session.binder.scaffold_pdb = f.read()
            session.binder.scaffold_file_path = scaffold_pdb
            session.update_stage_status(WorkflowStage.BINDER_SCAFFOLD_DESIGN, StageStatus.COMPLETED)
        
        # Load MPNN sequences (Step 3)
        mpnn_fasta = os.path.join(folder_path, f"{project_name}_Protein_MPNN_prediction.fa")
        if os.path.exists(mpnn_fasta):
            with open(mpnn_fasta, 'r') as f:
                fasta_content = f.read()
            session.binder.mpnn_fasta_content = fasta_content
            
            # Parse sequences
            sequences = []
            scores = []
            current_seq = None
            current_score = None
            
            for line in fasta_content.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_seq and current_score:
                        sequences.append(current_seq)
                        scores.append(float(current_score))
                    # Parse header: >seq_1, score=2.5862, global_score=2.5862
                    if 'score=' in line:
                        score_part = line.split('score=')[1].split(',')[0]
                        current_score = score_part
                    current_seq = ""
                elif line:
                    current_seq = line
            
            # Add last sequence
            if current_seq and current_score:
                sequences.append(current_seq)
                scores.append(float(current_score))
            
            session.binder.mpnn_sequences = sequences
            session.binder.mpnn_scores = scores
            session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.COMPLETED)
        
        # Load complex structure (Step 4)
        complex_pdb = os.path.join(folder_path, f"{project_name}_pdb_1_MULTIMER.pdb")
        if os.path.exists(complex_pdb):
            with open(complex_pdb, 'r') as f:
                session.complex.complex_pdb = f.read()
            session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.COMPLETED)
        
        # Load pLDDT scores (Step 5)
        plddt_file = os.path.join(folder_path, "pLDDT_scores.txt")
        if os.path.exists(plddt_file):
            with open(plddt_file, 'r') as f:
                content = f.read().strip()
                try:
                    session.complex.plddt_score = float(content)
                except ValueError:
                    # File might contain filename and score: "filename.pdb  score"
                    parts = content.split()
                    for part in parts:
                        try:
                            session.complex.plddt_score = float(part.strip())
                            break
                        except ValueError:
                            continue
        
        # Update project name
        session.project_name = project_name
        
        # Mark TARGET_INPUT as completed since we loaded a project
        session.update_stage_status(WorkflowStage.TARGET_INPUT, StageStatus.COMPLETED)
        
        # If we have complex, analyze the binding interface
        if session.complex.complex_pdb and session.target.pdb_content:
            try:
                # Parse structures
                target_atoms = parse_pdb_content(session.target.pdb_content)
                binder_pdb = session.binder.scaffold_pdb or session.binder.pdb_content
                
                if binder_pdb:
                    binder_atoms = parse_pdb_content(binder_pdb)
                    
                    # Find interface residues
                    target_interface = find_interface_residues(target_atoms, binder_atoms)
                    binder_interface = find_interface_residues(binder_atoms, target_atoms)
                    
                    # Assess binding quality
                    quality = assess_binding_quality(
                        target_atoms, binder_atoms,
                        target_interface, binder_interface
                    )
                    
                    # Update session with analysis
                    session.complex.target_interface_residues = target_interface
                    session.complex.binder_interface_residues = binder_interface
                    session.complex.num_contacts = quality['num_contacts']
                    session.complex.avg_distance = quality['avg_distance']
                    session.complex.min_distance = quality['min_distance']
                    session.complex.quality_score = quality['score']
                    session.complex.quality_grade = quality['grade']
                    session.complex.feedback = quality['feedback']
            except Exception as e:
                st.warning(f"⚠️ Could not analyze binding interface: {str(e)}")
        
        # Mark RESULTS as completed if we have a complex structure (regardless of analysis)
        if session.complex.complex_pdb:
            session.update_stage_status(WorkflowStage.RESULTS, StageStatus.COMPLETED)
        
        # Advance to appropriate stage based on what was loaded
        if session.complex.complex_pdb:
            session.advance_to_stage(WorkflowStage.RESULTS)
        elif session.binder.mpnn_sequences:
            session.advance_to_stage(WorkflowStage.COMPLEX_PREDICTION)
        elif session.binder.scaffold_pdb:
            session.advance_to_stage(WorkflowStage.BINDER_SEQUENCE_DESIGN)
        elif session.target.structure_predicted:
            session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
        
        return True
        
    except Exception as e:
        st.error(f"Error loading pipeline results: {str(e)}")
        import traceback
        st.code(traceback.format_exc())
        return False


def render_target_input_stage():
    """Stage 1: Target Protein Input"""
    st.header("1️⃣ Target Protein Input")
    
    # Page description
    st.markdown("""
    <div style="background-color: #f0f8ff; padding: 15px; border-radius: 8px; border-left: 4px solid #76B900; margin-bottom: 20px;">
        <h4 style="margin: 0 0 10px 0; color: #1A1A1A;">📌 What is a Target Protein?</h4>
        <p style="margin: 0; color: #333;">
            The <b>target protein</b> is the molecule you want to design a binder for. This could be:
        </p>
        <ul style="margin: 10px 0; color: #333;">
            <li><b>A disease-related protein</b> - e.g., a viral spike protein you want to neutralize</li>
            <li><b>A receptor</b> - e.g., a cell surface receptor you want to activate or block</li>
            <li><b>An enzyme</b> - e.g., an enzyme you want to inhibit</li>
            <li><b>An antibody target</b> - e.g., an antigen for therapeutic antibody development</li>
        </ul>
        <p style="margin: 0; color: #666; font-size: 0.9em;">
            💡 <i>Example: To design a binder for insulin, you would input the insulin sequence as your target.</i>
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    session = st.session_state.workflow_session
    target = session.target
    
    # Project mode selection
    st.subheader("Project Mode")
    project_mode = st.radio(
        "Choose how to start",
        ["Create New Project", "Load Previous Results"],
        horizontal=True,
        key="project_mode"
    )
    
    if project_mode == "Load Previous Results":
        st.markdown("---")
        st.markdown("### Load from existing pipeline output")
        
        # Scan for available output folders
        import glob
        workflow_outputs = glob.glob("workflow/*_output")
        root_outputs = glob.glob("*_output")
        all_outputs = workflow_outputs + root_outputs
        
        if all_outputs:
            # Create display names with metadata
            output_options = {}
            for folder in all_outputs:
                folder_name = os.path.basename(folder)
                if '_AF2_output' in folder_name:
                    project = folder_name.replace('_AF2_output', '')
                    display = f"{project} (AlphaFold2)"
                elif '_OF3_output' in folder_name:
                    project = folder_name.replace('_OF3_output', '')
                    display = f"{project} (OpenFold3)"
                elif '_AF3_output' in folder_name:
                    project = folder_name.replace('_AF3_output', '')
                    display = f"{project} (AlphaFold3)"
                else:
                    project = folder_name.replace('_output', '')
                    display = f"{project}"
                output_options[display] = folder
            
            selected_display = st.selectbox(
                "Select Previous Project",
                options=list(output_options.keys()),
                help="Choose from available pipeline output folders"
            )
            
            selected_folder = output_options[selected_display]
            
            # Show folder details
            col1, col2 = st.columns(2)
            with col1:
                st.info(f"Folder: `{selected_folder}`")
            with col2:
                # Check what stages are available
                folder_name = os.path.basename(selected_folder)
                project = folder_name.replace('_AF2_output', '').replace('_OF3_output', '').replace('_output', '')
                has_target = os.path.exists(os.path.join(selected_folder, f"{project}_first_structure.pdb"))
                has_scaffold = os.path.exists(os.path.join(selected_folder, f"{project}_RFD_prediction.pdb"))
                has_sequences = os.path.exists(os.path.join(selected_folder, f"{project}_Protein_MPNN_prediction.fa"))
                has_complex = os.path.exists(os.path.join(selected_folder, f"{project}_pdb_1_MULTIMER.pdb"))
                
                stages_text = []
                if has_target: stages_text.append("Target")
                if has_scaffold: stages_text.append("Scaffold")
                if has_sequences: stages_text.append("Sequences")
                if has_complex: stages_text.append("Complex")
                
                st.success(f"✅ Stages: {', '.join(stages_text)}")
            
            if st.button("Load This Project", type="primary", use_container_width=True):
                with st.spinner("Loading pipeline results..."):
                    load_success = load_pipeline_results(session, selected_folder)
                    if load_success:
                        st.success("✅ Successfully loaded pipeline results!")
                        st.balloons()
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error("❌ Failed to load results. Check folder structure.")
        else:
            st.warning("⚠️ No previous pipeline output folders found in workspace")
            st.info("Create a new project to get started!")
        
        return  # Exit early when in load mode
    
    # Create New Project mode
    st.markdown("---")
    st.markdown("### Define the protein you want to design a binder for")
    
    # Input type selection
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Map input_type to display options
        input_type_map = {
            "sequence": "Sequence",
            "pdb_file": "PDB File",
            "pdb_id": "PDB ID"
        }
        reverse_map = {v: k for k, v in input_type_map.items()}
        
        # Get current index
        current_display = input_type_map.get(target.input_type, "Sequence")
        current_index = ["Sequence", "PDB File", "PDB ID"].index(current_display)
        
        input_type = st.radio(
            "Input Type",
            ["Sequence", "PDB File", "PDB ID"],
            horizontal=True,
            index=current_index
        )
        target.input_type = reverse_map[input_type]
    
    # Input fields based on type
    if input_type == "Sequence":
        sequence = st.text_area(
            "Target Protein Sequence",
            value=target.sequence or "",
            placeholder="Enter amino acid sequence (e.g., MKTAYIAKQRQISFVKSHF...)",
            height=150,
            help="Enter the amino acid sequence of your target protein"
        )
        
        if sequence:
            is_valid, result = validate_protein_sequence(sequence)
            if is_valid:
                target.sequence = result
                st.success(f"✅ Valid sequence: {len(result)} amino acids")
            else:
                st.error(f"❌ {result}")
    
    elif input_type == "PDB File":
        uploaded_file = st.file_uploader("Upload PDB File", type=['pdb'])
        if uploaded_file:
            pdb_content = uploaded_file.read().decode('utf-8')
            target.pdb_content = pdb_content
            target.structure_predicted = True
            target.model_used = "Uploaded PDB"
            st.success("✅ PDB file loaded successfully")
            
            # Extract sequence from PDB
            atoms = parse_pdb_content(pdb_content)
            if atoms:
                st.info(f"Structure contains {len(atoms)} atoms")
    
    elif input_type == "PDB ID":
        pdb_id = st.text_input(
            "PDB ID",
            value=target.pdb_id or "",
            placeholder="e.g., 7BV2",
            max_chars=4,
            help="Enter a 4-character PDB ID from RCSB PDB"
        )
        
        if pdb_id and len(pdb_id) == 4:
            target.pdb_id = pdb_id.upper()
            if st.button("Fetch from PDB"):
                with st.spinner("Fetching structure..."):
                    # Implement PDB fetch logic here
                    st.info("PDB fetching not yet implemented")
    
    # Binding site specification (optional)
    with st.expander("Advanced: Specify Binding Site (Optional)"):
        st.markdown("Specify residues that should be part of the binding interface")
        binding_site = st.text_input(
            "Binding Site Residues",
            placeholder="e.g., 10-20, 45, 67-72",
            help="Comma-separated residue numbers or ranges"
        )
        
        if binding_site:
            # Parse binding site input
            try:
                residues = []
                for part in binding_site.split(','):
                    part = part.strip()
                    if '-' in part:
                        start, end = map(int, part.split('-'))
                        residues.extend(range(start, end + 1))
                    else:
                        residues.append(int(part))
                target.binding_site_residues = sorted(set(residues))
                st.success(f"✅ Binding site: {len(target.binding_site_residues)} residues")
            except ValueError:
                st.error("❌ Invalid format. Use: 10-20, 45, 67-72")
    
    # Navigation
    col1, col2, col3 = st.columns([1, 1, 2])
    
    with col2:
        if st.button("Next: Predict Structure →", type="primary", use_container_width=True):
            can_advance, message = session.can_advance_to(WorkflowStage.TARGET_PREDICTION)
            if can_advance:
                session.advance_to_stage(WorkflowStage.TARGET_PREDICTION)
                session.update_stage_status(WorkflowStage.TARGET_INPUT, StageStatus.COMPLETED)
                st.rerun()
            else:
                st.error(f"❌ {message}")


def render_target_prediction_stage():
    """Stage 2: Target Structure Prediction"""
    st.header("2️⃣ Target Structure Prediction")
    
    # Page description
    st.markdown("""
    <div style="background-color: #f0f8ff; padding: 15px; border-radius: 8px; border-left: 4px solid #76B900; margin-bottom: 20px;">
        <p style="margin: 0; color: #333;">
            Predict the 3D structure of your target protein using state-of-the-art AI models. 
            The structure is essential for designing a binder that fits properly.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    session = st.session_state.workflow_session
    target = session.target
    
    # Get pipeline (lazy loading) - use GenerativePipeline for better model handling
    def get_pipeline():
        if 'pipeline' not in st.session_state:
            st.session_state.pipeline = GenerativePipeline(
                session=session,
                api_key=st.session_state.api_key,
                output_dir=Path(f"{session.project_name}_output")
            )
        return st.session_state.pipeline
    
    # Show target summary
    st.subheader("Target Summary")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if target.sequence:
            st.metric("Sequence Length", f"{len(target.sequence)} AA")
    with col2:
        st.metric("Input Type", target.input_type.replace("_", " ").title())
    with col3:
        st.metric("Status", "✅ Ready" if target.sequence or target.pdb_content else "⏳ Pending")
    
    # Skip if structure already available
    if target.pdb_content:
        st.info("✅ Structure already available (uploaded PDB)")
        session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
    else:
        # Model selection with enhanced options
        st.subheader("Structure Prediction Model")
        
        # Define models with detailed info
        model_info = {
            "AlphaFold3 (Self-hosted)": {
                "id": "AF3",
                "speed": "🚀 1-5 min",
                "accuracy": "⭐⭐⭐⭐⭐",
                "desc": "Latest DeepMind model - hosted at brevlab.com",
                "supports_algorithm": False,
                "supports_diffusion": True
            },
            "AlphaFold2": {
                "id": "AF2",
                "speed": "🐌 5-10 min",
                "accuracy": "⭐⭐⭐⭐⭐",
                "desc": "High accuracy, multiple structure predictions (NVIDIA API)",
                "supports_algorithm": True,
                "supports_diffusion": False
            },
            "OpenFold3": {
                "id": "OF3",
                "speed": "⚡ 2-5 min",
                "accuracy": "⭐⭐⭐⭐",
                "desc": "Fast and accurate, good for prototyping (NVIDIA API)",
                "supports_algorithm": False,
                "supports_diffusion": False
            }
        }
        
        # Model selection
        col1, col2 = st.columns([2, 1])
        
        with col1:
            selected_model_name = st.selectbox(
                "Select Prediction Model",
                options=list(model_info.keys()),
                index=0,  # Default to AlphaFold3
                help="AlphaFold3 is the latest model (self-hosted), AlphaFold2/OpenFold3 use NVIDIA API"
            )
        
        with col2:
            info = model_info[selected_model_name]
            st.markdown(f"""
            **Speed:** {info['speed']}  
            **Accuracy:** {info['accuracy']}
            """)
        
        st.caption(f"💡 {info['desc']}")
        
        # Model-specific options
        algorithm = "mmseqs2"  # default for AF2
        num_diffusion_samples = 1  # default for AF3
        model_seeds = [42]  # default for AF3
        
        if info.get("supports_diffusion"):
            with st.expander("⚙️ AlphaFold3 Options", expanded=False):
                st.markdown("### Diffusion Sampling Parameters")
                num_diffusion_samples = st.slider(
                    "Number of Diffusion Samples",
                    min_value=1,
                    max_value=5,
                    value=1,
                    help="More samples = better coverage but slower"
                )
                seed_input = st.number_input(
                    "Random Seed",
                    min_value=1,
                    max_value=9999,
                    value=42,
                    help="Set seed for reproducible results"
                )
                model_seeds = [seed_input]
                
                st.info("🔬 AlphaFold3 uses diffusion-based structure prediction for improved accuracy on challenging targets.")
        
        if info.get("supports_algorithm"):
            with st.expander("⚙️ Advanced Model Options", expanded=False):
                st.markdown("### AlphaFold2 MSA Algorithm")
                algorithm = st.radio(
                    "MSA Algorithm",
                    ["mmseqs2", "jackhmmer"],
                    index=0,
                    horizontal=True,
                    help="mmseqs2 is faster, jackhmmer is more sensitive for distant homologs"
                )
                
                st.info("""
                **mmseqs2** (Recommended): Fast MSA search, works well for most proteins  
                **jackhmmer**: Slower but more thorough, better for sequences with few homologs
                """)
        
        # Show processing time warning for long sequences
        if target.sequence and len(target.sequence) > 300:
            st.warning(f"⚠️ Your sequence has {len(target.sequence)} residues. {selected_model_name} may take longer to process.")
        
        # Predict button
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            predict_clicked = st.button(
                f"🔬 Predict with {selected_model_name}", 
                type="primary",
                use_container_width=True
            )
        
        if predict_clicked:
            if not st.session_state.api_key and not st.session_state.demo_mode:
                st.error("❌ Please provide an API key or enable Demo Mode")
                return
            
            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.IN_PROGRESS)
            
            # Show processing info
            if selected_model_name == "AlphaFold2":
                st.info("⏳ AlphaFold2 typically takes 5-10 minutes. The model will generate 5 structure predictions ranked by confidence.")
            
            if info["id"] == "AF3":
                st.info("🔬 Connecting to AlphaFold3 server at brevlab.com...")
            
            with st.spinner(f"Predicting structure with {selected_model_name}..."):
                try:
                    if st.session_state.demo_mode:
                        # Demo mode
                        time.sleep(2)
                        from frontend.app_v2 import generate_mock_pdb
                        pdb_content = generate_mock_pdb(target.sequence)
                        target.pdb_content = pdb_content
                        target.structure_predicted = True
                        target.model_used = selected_model_name
                    else:
                        # Use GenerativePipeline for real prediction
                        pipeline = get_pipeline()
                        
                        # Run target prediction via pipeline
                        # Pass different params based on model type
                        if info["id"] == "AF3":
                            success, msg = pipeline.run_target_prediction(
                                model=info["id"],
                                num_diffusion_samples=num_diffusion_samples,
                                model_seeds=model_seeds
                            )
                        else:
                            success, msg = pipeline.run_target_prediction(
                                model=info["id"],
                                algorithm=algorithm
                            )
                        
                        if not success:
                            st.error(f"❌ Prediction failed: {msg}")
                            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
                            return
                    
                    session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
                    st.success(f"🎉 Structure prediction completed with {selected_model_name}!")
                    st.balloons()
                    time.sleep(1)
                    st.rerun()
                    
                except Exception as e:
                    st.error(f"❌ Error: {str(e)}")
                    import traceback
                    st.code(traceback.format_exc())
                    session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
    
    # Show structure if available
    if target.pdb_content:
        st.subheader("Target Structure")
        
        # Calculate average pLDDT from PDB if not already set
        avg_plddt = target.confidence_avg
        if not avg_plddt and target.pdb_content:
            b_factors = []
            for line in target.pdb_content.split('\\n'):
                if line.startswith('ATOM'):
                    try:
                        b_factor = float(line[60:66].strip())
                        b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
            if b_factors:
                # Check if B-factors look like pLDDT scores
                min_b, max_b = min(b_factors), max(b_factors)
                if 0 <= min_b <= 100 and 0 <= max_b <= 100 and (max_b - min_b) > 5:
                    avg_plddt = sum(b_factors) / len(b_factors)
        
        # Show model info
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Model Used", target.model_used or "Unknown")
        with col2:
            if avg_plddt:
                st.metric("Avg pLDDT", f"{avg_plddt:.1f}")
            else:
                st.metric("Structure", "✅ Predicted")
        with col3:
            if target.all_structures_pdb:
                st.metric("Predictions", "5 models (AF2)")
            else:
                st.metric("Predictions", "1 model")
        
        # If AlphaFold2 returned multiple structures, show selector
        if target.all_structures_pdb and "MODEL" in target.all_structures_pdb:
            with st.expander("📊 View All 5 AlphaFold2 Predictions", expanded=False):
                st.info("AlphaFold2 generates 5 independent structure predictions ranked by confidence. The top-ranked structure is shown by default.")
                
                # Parse multiple models from all_structures_pdb
                models = target.all_structures_pdb.split("ENDMDL")
                num_models = len([m for m in models if m.strip()])
                
                selected_model = st.selectbox(
                    "Select Model",
                    options=list(range(1, min(6, num_models + 1))),
                    format_func=lambda x: f"Model {x} {'(Top Ranked)' if x == 1 else ''}",
                    help="View different AlphaFold2 structure predictions"
                )
                
                st.caption(f"Showing model {selected_model} of {num_models}")
        
        # 3D Visualization with pLDDT coloring and vertical legend bar
        try:
            has_plddt = check_has_plddt_scores(target.pdb_content)
            html_content = create_3d_visualization(target.pdb_content, color_by_plddt=has_plddt)
            st.components.v1.html(html_content, height=620)
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
        
        # Download buttons
        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                "📥 Download Target PDB",
                data=target.pdb_content,
                file_name=f"target_{session.project_name}.pdb",
                mime="chemical/x-pdb"
            )
        with col2:
            if target.all_structures_pdb:
                st.download_button(
                    "📥 Download All 5 Models",
                    data=target.all_structures_pdb,
                    file_name=f"target_{session.project_name}_all_models.pdb",
                    mime="chemical/x-pdb"
                )
    
    # Navigation
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back", use_container_width=True):
            session.advance_to_stage(WorkflowStage.TARGET_INPUT)
            st.rerun()
    
    with col2:
        if st.button("Next: Design Binder →", type="primary", use_container_width=True):
            if target.pdb_content:
                session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
                st.rerun()
            else:
                st.error("❌ Please predict target structure first")


def render_binder_design_stage():
    """DEPRECATED - Replaced by render_binder_scaffold_stage"""
    render_binder_scaffold_stage()


def render_binder_scaffold_stage():
    """Stage 3: Binder Scaffold Design (RFDiffusion)"""
    st.header("3️⃣ Binder Scaffold Design")
    
    # Page description
    st.markdown("""
    <div style="background-color: #f0f8ff; padding: 15px; border-radius: 8px; border-left: 4px solid #76B900; margin-bottom: 20px;">
        <p style="margin: 0; color: #333;">
            <b>RFDiffusion</b> uses AI diffusion models to generate a 3D backbone scaffold for your binder protein. 
            This scaffold defines the overall shape and structure that will interact with your target.
        </p>
        <p style="margin: 10px 0 0 0; color: #666; font-size: 0.9em;">
            💡 <i>The scaffold is like an architectural blueprint - it defines the shape but not the amino acid sequence.</i>
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    session = st.session_state.workflow_session
    binder = session.binder
    
    # Initialize pipeline only when needed (lazy loading)
    def get_pipeline():
        if 'pipeline' not in st.session_state:
            st.session_state.pipeline = GenerativePipeline(
                session=session,
                api_key=st.session_state.api_key,
                output_dir=Path(f"{session.project_name}_output")
            )
        return st.session_state.pipeline
    
    # Show target summary
    st.subheader("Target Protein")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Length", f"{len(session.target.sequence)} AA")
    with col2:
        st.metric("Model Used", session.target.model_used or "N/A")
    with col3:
        st.metric("Structure", "✅ Available" if session.target.pdb_content else "❌ Missing")
    
    # Extract available residues from target PDB for hotspot guidance
    available_residues = {}
    if session.target.pdb_content:
        from workflow.generative_pipeline import extract_residues_from_pdb
        available_residues = extract_residues_from_pdb(session.target.pdb_content)
    
    # RFDiffusion parameters
    st.subheader("RFDiffusion Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Show available residue ranges to help with contigs
        default_contigs = "A1-25/0 70-100"  # Fallback default
        if available_residues:
            residue_info = ", ".join([f"Chain {chain}: {min(nums)}-{max(nums)}" for chain, nums in available_residues.items()])
            st.caption(f"📊 Available residues: {residue_info}")
            
            # Generate a valid default contigs based on actual PDB residues
            for chain, nums in available_residues.items():
                if len(nums) >= 10:
                    min_res = min(nums)
                    # Use first 25 residues or fewer if chain is shorter
                    end_res = min(min_res + 24, max(nums))
                    default_contigs = f"{chain}{min_res}-{end_res}/0 70-100"
                    break
        
        # Use stored value if available, otherwise use computed default
        stored_contigs = binder.rfdiffusion_params.get("contigs", "")
        contigs = st.text_input(
            "Contigs Specification",
            value=stored_contigs if stored_contigs else default_contigs,
            help="Format: A1-25/0 70-100 means keep residues 1-25 from target, add 70-100 new residues for binder"
        )
        
        # Validate contigs in real-time
        if available_residues and contigs:
            from workflow.generative_pipeline import validate_and_fix_contigs
            fixed, warnings = validate_and_fix_contigs(session.target.pdb_content, contigs)
            if warnings:
                st.warning(f"⚠️ {'; '.join(warnings)}")
                if fixed != contigs:
                    st.info(f"💡 Suggested fix: `{fixed}`")
        
        diffusion_steps = st.slider(
            "Diffusion Steps",
            min_value=5,
            max_value=50,
            value=15,
            help="More steps = better quality but slower"
        )
    
    with col2:
        # Show example hotspots based on actual available residues
        example_hotspots = ""
        if available_residues:
            for chain, nums in available_residues.items():
                if len(nums) >= 4:
                    # Pick some residues from the middle of the chain
                    mid = len(nums) // 2
                    example_hotspots = f"{chain}{nums[mid]}, {chain}{nums[mid+1]}, {chain}{nums[mid+2]}"
                    break
        
        hotspot_input = st.text_input(
            "Binding Hotspot Residues (optional)",
            value=", ".join(binder.rfdiffusion_params.get("hotspot_res", [])) if binder.rfdiffusion_params.get("hotspot_res") else "",
            placeholder=f"e.g., {example_hotspots}" if example_hotspots else "e.g., A14, A15, A17",
            help="Specify target residues that must be in the binding interface. Format: ChainResidue (e.g., A14)"
        )
        
        hotspots = [h.strip() for h in hotspot_input.split(",") if h.strip()] if hotspot_input else None
        
        # Validate hotspots if entered
        if hotspots and session.target.pdb_content:
            from workflow.generative_pipeline import validate_hotspot_residues
            valid, invalid = validate_hotspot_residues(session.target.pdb_content, hotspots)
            if invalid:
                st.warning(f"⚠️ Invalid residues: {', '.join(invalid)}")
            if valid:
                st.success(f"✅ Valid hotspots: {', '.join(valid)}")
    
    # Check if stage is completed (view-only mode)
    stage_status = session.stage_statuses.get(WorkflowStage.BINDER_SCAFFOLD_DESIGN.value)
    is_completed = stage_status == StageStatus.COMPLETED.value
    
    # Display existing scaffold if available
    if binder.scaffold_pdb:
        st.success("✅ Binder scaffold already generated!")
        
        with st.expander("View Scaffold Structure", expanded=True):
            # RFDiffusion scaffold doesn't have pLDDT - use single color
            viz_html = create_3d_visualization(binder.scaffold_pdb, color_by_plddt=False)
            st.components.v1.html(viz_html, height=600, scrolling=False)
        
        if is_completed:
            st.info("💡 This stage is completed. Use the navigation buttons below to proceed or click 'Regenerate' to create a new scaffold.")
    
    # Generate/Regenerate button - only show if not in view-only mode OR explicitly want to regenerate
    if not is_completed:
        col1, col2, col3 = st.columns([1, 2, 1])
        
        with col2:
            if st.button("Generate Binder Scaffold", type="primary", use_container_width=True):
                pipeline = get_pipeline()
                
                # Show debug info
                st.info("Debug: Running RFDiffusion... Check your terminal for detailed logs!")
                
                with st.spinner("Running RFDiffusion... This may take a few minutes..."):
                    success, msg = pipeline.run_scaffold_design(
                        contigs=contigs,
                        hotspot_res=hotspots,
                        diffusion_steps=diffusion_steps
                    )
                    
                    if success:
                        st.success(f"✅ {msg}")
                        st.balloons()
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error(f"❌ {msg}")
    else:
        # Show regenerate option for completed stage
        with st.expander("Regenerate with Different Parameters"):
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                if st.button("Regenerate Scaffold", type="secondary", use_container_width=True):
                    # Reset status and allow regeneration
                    session.update_stage_status(WorkflowStage.BINDER_SCAFFOLD_DESIGN, StageStatus.IN_PROGRESS)
                    st.rerun()
    
    # Navigation
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back to Target", use_container_width=True):
            session.advance_to_stage(WorkflowStage.TARGET_PREDICTION)
            st.rerun()
    
    with col3:
        if st.button("Next: Sequence Design →", type="primary", use_container_width=True, disabled=not binder.scaffold_pdb):
            session.advance_to_stage(WorkflowStage.BINDER_SEQUENCE_DESIGN)
            st.rerun()


def render_binder_sequence_stage():
    """Stage 4: Sequence Design (ProteinMPNN)"""
    st.header("4️⃣ Sequence Design")
    
    # Page description
    st.markdown("""
    <div style="background-color: #f0f8ff; padding: 15px; border-radius: 8px; border-left: 4px solid #76B900; margin-bottom: 20px;">
        <p style="margin: 0; color: #333;">
            <b>ProteinMPNN</b> designs amino acid sequences that will fold into your scaffold structure. 
            It generates multiple sequence variants optimized for stability and binding.
        </p>
        <p style="margin: 10px 0 0 0; color: #666; font-size: 0.9em;">
            💡 <i>Think of it as filling in the blueprint with specific building materials (amino acids).</i>
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    session = st.session_state.workflow_session
    binder = session.binder
    
    # Check if stage is completed (view-only mode)
    stage_status = session.stage_statuses.get(WorkflowStage.BINDER_SEQUENCE_DESIGN.value)
    is_completed = stage_status == StageStatus.COMPLETED.value
    
    # Get pipeline (lazy loading)
    def get_pipeline():
        if 'pipeline' not in st.session_state:
            st.session_state.pipeline = GenerativePipeline(
                session=session,
                api_key=st.session_state.api_key,
                output_dir=Path(f"{session.project_name}_output")
            )
        return st.session_state.pipeline
    
    # Show scaffold summary
    st.subheader("Binder Scaffold")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Scaffold", "✅ Available" if binder.scaffold_pdb else "❌ Missing")
    with col2:
        # Show the actual model used for scaffold design
        scaffold_method = binder.rfdiffusion_params.get("model", "RFDiffusion") if binder.rfdiffusion_params else "RFDiffusion"
        st.metric("Design Method", scaffold_method)
    with col3:
        num_existing = len(binder.mpnn_sequences)
        st.metric("Sequences", f"{num_existing} designed" if num_existing > 0 else "Not generated")
    
    # ProteinMPNN parameters
    st.subheader("ProteinMPNN Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        num_sequences = st.slider(
            "Number of Sequences",
            min_value=1,
            max_value=20,
            value=10,
            help="How many sequence variants to generate"
        )
    
    with col2:
        sampling_temp = st.slider(
            "Sampling Temperature",
            min_value=0.01,
            max_value=0.5,
            value=0.1,
            step=0.01,
            help="Lower = more conservative, Higher = more diverse"
        )
    
    # Display existing sequences
    if binder.mpnn_sequences:
        # Filter out sequences that are entirely glycine linkers
        import re
        valid_sequences = []
        for i, seq in enumerate(binder.mpnn_sequences):
            # Remove runs of 5+ glycines
            display_seq = re.sub(r'G{5,}', '', seq)
            # Skip if empty or only has a few amino acids left (likely just linker)
            if len(display_seq.strip()) >= 10:
                score = binder.mpnn_scores[i] if i < len(binder.mpnn_scores) else None
                valid_sequences.append((i, seq, display_seq, score))
        
        st.success(f"✅ {len(valid_sequences)} valid sequences (filtered from {len(binder.mpnn_sequences)} total)")
        
        with st.expander("📋 View Designed Sequences", expanded=True):
            for orig_idx, seq, display_seq, score in valid_sequences:
                selected = "[*]" if orig_idx == binder.selected_sequence_idx else "[ ]"
                
                col1, col2, col3 = st.columns([1, 8, 2])
                with col1:
                    st.markdown(f"**{selected} #{orig_idx+1}**")
                with col2:
                    # Show cleaned sequence (without linkers)
                    if display_seq != seq:
                        st.caption("(Linker removed)")
                    st.code(display_seq[:80] + "..." if len(display_seq) > 80 else display_seq, language="text")
                with col3:
                    if score:
                        st.metric("Score", f"{score:.3f}")
                    if st.button(f"Select", key=f"select_{orig_idx}"):
                        binder.selected_sequence_idx = orig_idx
                        binder.sequence = seq  # Store original sequence with linker
                        st.success(f"Selected sequence #{orig_idx+1}")
                        st.rerun()
        
        if is_completed:
            st.info("💡 This stage is completed. Use the navigation buttons below to proceed or click 'Regenerate' to design new sequences.")
        else:
            st.info(f"💡 Currently selected: Sequence #{binder.selected_sequence_idx + 1}")
    
    # Generate/Regenerate button
    if not is_completed:
        col1, col2, col3 = st.columns([1, 2, 1])
        
        with col2:
            if st.button("Design Sequences", type="primary", use_container_width=True):
                pipeline = get_pipeline()
                
                # Show debug info
                st.info("Debug: Running ProteinMPNN... Check your terminal for detailed logs!")
                
                with st.spinner("Running ProteinMPNN... Designing sequences..."):
                    success, msg = pipeline.run_sequence_design(
                        num_sequences=num_sequences,
                        sampling_temp=sampling_temp
                    )
                    
                    if success:
                        st.success(f"✅ {msg}")
                        st.balloons()
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error(f"❌ {msg}")
    else:
        # Show regenerate option for completed stage
        with st.expander("Regenerate with Different Parameters"):
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                if st.button("Regenerate Sequences", type="secondary", use_container_width=True):
                    # Reset status and allow regeneration
                    session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.IN_PROGRESS)
                    st.rerun()
    
    # Navigation
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back to Scaffold", use_container_width=True):
            session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
            st.rerun()
    
    with col3:
        if st.button("Next: Predict Complex →", type="primary", use_container_width=True, disabled=not binder.mpnn_sequences):
            session.advance_to_stage(WorkflowStage.COMPLEX_PREDICTION)
            st.rerun()


def render_binder_prediction_stage():
    """DEPRECATED - Replaced by render_complex_prediction_stage"""
    render_complex_prediction_stage()


def render_complex_prediction_stage():
    """Stage 5: Complex Prediction (AlphaFold-Multimer)"""
    st.header("5️⃣ Complex Prediction")
    
    # Page description
    st.markdown("""
    <div style="background-color: #f0f8ff; padding: 15px; border-radius: 8px; border-left: 4px solid #76B900; margin-bottom: 20px;">
        <p style="margin: 0; color: #333;">
            <b>AlphaFold-Multimer</b> predicts how your designed binder will interact with the target protein. 
            This shows the 3D structure of the binder-target complex and validates the binding interaction.
        </p>
        <p style="margin: 10px 0 0 0; color: #666; font-size: 0.9em;">
            💡 <i>A successful prediction shows the binder positioned at the target's binding site with good confidence scores.</i>
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    session = st.session_state.workflow_session
    binder = session.binder
    complex_data = session.complex
    
    # Check if stage is completed (view-only mode)
    stage_status = session.stage_statuses.get(WorkflowStage.COMPLEX_PREDICTION.value)
    is_completed = stage_status == StageStatus.COMPLETED.value
    
    # Get pipeline (lazy loading)
    def get_pipeline():
        if 'pipeline' not in st.session_state:
            st.session_state.pipeline = GenerativePipeline(
                session=session,
                api_key=st.session_state.api_key,
                output_dir=Path(f"{session.project_name}_output")
            )
        return st.session_state.pipeline
    
    # Show summary
    st.subheader("Complex Components")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        target_seq = session.target.sequence or ""
        st.metric("Target", f"{len(target_seq)} AA")
    with col2:
        # Get binder sequence from either sequence field or mpnn_sequences
        binder_seq = binder.sequence or (binder.mpnn_sequences[binder.selected_sequence_idx] if binder.mpnn_sequences else "")
        st.metric("Binder", f"{len(binder_seq)} AA")
    with col3:
        st.metric("Selected Seq", f"#{binder.selected_sequence_idx + 1}")
    with col4:
        st.metric("Complex", "✅ Predicted" if complex_data.complex_pdb else "Not predicted")
    
    # Complex Prediction Parameters
    st.subheader("Complex Prediction Parameters")
    
    # Model selection - AF2-Multimer vs AF3
    complex_model_options = {
        "alphafold2_multimer": {
            "name": "🔬 AlphaFold2-Multimer",
            "desc": "NVIDIA BioNeMo API - optimized for protein complexes"
        },
        "alphafold3": {
            "name": "🧬 AlphaFold3",
            "desc": "Latest DeepMind model - hosted at brevlab.com"
        }
    }
    
    selected_complex_model = st.radio(
        "Select Model",
        options=list(complex_model_options.keys()),
        format_func=lambda x: complex_model_options[x]["name"],
        horizontal=True,
        help="Choose which model to use for complex structure prediction"
    )
    
    st.caption(f"*{complex_model_options[selected_complex_model]['desc']}*")
    
    col1, col2 = st.columns(2)
    
    with col1:
        num_candidates = st.slider(
            "Number of Candidates to Evaluate",
            min_value=1,
            max_value=min(10, len(binder.mpnn_sequences)),
            value=min(3, len(binder.mpnn_sequences)),
            help="Predict complexes for top N designed sequences"
        )
    
    with col2:
        if selected_complex_model == "alphafold2_multimer":
            selected_models = st.multiselect(
                "AlphaFold2 Models",
                options=[1, 2, 3, 4, 5],
                default=[1],
                help="Which AF-Multimer models to use (more = slower but more accurate)"
            )
        else:
            st.info("AlphaFold3 uses a single unified model")
            selected_models = [1]  # Placeholder for AF3
    
    # Display existing complex
    if complex_data.complex_pdb:
        st.success(f"✅ Complex predicted! pLDDT Score: {complex_data.plddt_score:.2f}")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            with st.expander("View Complex Structure", expanded=True):
                viz_html = create_3d_visualization(complex_data.complex_pdb)
                st.components.v1.html(viz_html, height=600, scrolling=False)
        
        with col2:
            st.metric("pLDDT Score", f"{complex_data.plddt_score:.2f}")
            st.metric("Quality Grade", complex_data.quality_grade)
            st.metric("Method", complex_data.docking_method)
            
            if complex_data.plddt_score > 90:
                st.success("Excellent confidence!")
            elif complex_data.plddt_score > 70:
                st.info("Good confidence")
            else:
                st.warning("⚠️ Low confidence")
        
        # Show rankings if multiple candidates
        if complex_data.candidate_rankings:
            with st.expander(f"View All {len(complex_data.candidate_rankings)} Candidates"):
                for i, candidate in enumerate(complex_data.candidate_rankings):
                    col1, col2, col3, col4 = st.columns([1, 6, 2, 2])
                    with col1:
                        st.markdown(f"**#{i+1}**")
                    with col2:
                        st.code(candidate['sequence'][:60] + "...", language="text")
                    with col3:
                        st.metric("pLDDT", f"{candidate['plddt_score']:.2f}")
                    with col4:
                        st.markdown(f"*{candidate['quality_grade']}*")
        
        if is_completed:
            st.info("💡 This stage is completed. Use the navigation buttons below to proceed or click 'Regenerate' to create a new prediction.")
    
    # Predict button
    if not is_completed:
        col1, col2, col3 = st.columns([1, 2, 1])
        
        with col2:
            model_display_name = "AlphaFold3" if selected_complex_model == "alphafold3" else "AlphaFold2-Multimer"
            if num_candidates > 1:
                button_label = f"Predict {num_candidates} Complexes with {model_display_name}"
            else:
                button_label = f"Predict Complex with {model_display_name}"
            
            if st.button(button_label, type="primary", use_container_width=True):
                pipeline = get_pipeline()
                
                # Show debug info
                if selected_complex_model == "alphafold3":
                    st.info(f"🧬 Connecting to AlphaFold3 server at brevlab.com...")
                else:
                    st.info(f"🔬 Running AlphaFold2-Multimer for {num_candidates} candidate(s)...")
                
                with st.spinner(f"Running {model_display_name} for {num_candidates} candidate(s)... This will take several minutes..."):
                    if num_candidates > 1:
                        success, msg = pipeline.run_batch_complex_prediction(
                            num_candidates=num_candidates,
                            selected_models=selected_models,
                            model_type=selected_complex_model
                        )
                    else:
                        success, msg = pipeline.run_complex_prediction(
                            sequence_idx=binder.selected_sequence_idx,
                            selected_models=selected_models,
                            model_type=selected_complex_model
                        )
                    
                    if success:
                        st.success(f"✅ {msg}")
                        st.balloons()
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error(f"❌ {msg}")
    else:
        # Show regenerate option for completed stage
        with st.expander("Regenerate with Different Parameters"):
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                if st.button("Regenerate Complex", type="secondary", use_container_width=True):
                    # Reset status and allow regeneration
                    session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.IN_PROGRESS)
                    st.rerun()
    
    # Navigation
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back to Sequences", use_container_width=True):
            session.advance_to_stage(WorkflowStage.BINDER_SEQUENCE_DESIGN)
            st.rerun()
    
    with col3:
        if st.button("View Results →", type="primary", use_container_width=True, disabled=not complex_data.complex_pdb):
            session.advance_to_stage(WorkflowStage.RESULTS)
            st.rerun()


def render_complex_analysis_stage():
    """DEPRECATED - Merged into render_complex_prediction_stage"""
    render_complex_prediction_stage()


def render_binder_design_stage():
    """DEPRECATED - Replaced by render_binder_scaffold_stage"""
    render_binder_scaffold_stage()


def render_binder_prediction_stage():
    """Stage 4: Binder Structure Prediction"""
    st.header("4️⃣ Binder Structure Prediction")
    
    session = st.session_state.workflow_session
    binder = session.binder
    
    # Show binder summary
    st.subheader("Binder Summary")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Sequence Length", f"{len(binder.sequence)} AA")
    with col2:
        st.metric("Design Method", binder.design_method.title())
    with col3:
        st.metric("Status", "✅ Ready")
    
    # Model selection
    st.subheader("Structure Prediction")
    
    model_options = {name: model for name, model in PROTEIN_MODELS.items()}
    selected_model_name = st.selectbox(
        "Select Prediction Model",
        options=list(model_options.keys()),
        index=1
    )
    
    selected_model = model_options[selected_model_name]
    
    # Predict button
    if st.button("Predict Binder Structure", type="primary"):
        if not st.session_state.api_key and not st.session_state.demo_mode:
            st.error("❌ Please provide an API key or enable Demo Mode")
            return
        
        session.update_stage_status(WorkflowStage.BINDER_PREDICTION, StageStatus.IN_PROGRESS)
        
        with st.spinner(f"Predicting binder structure with {selected_model_name}..."):
            try:
                if st.session_state.demo_mode:
                    time.sleep(2)
                    from frontend.app_v2 import generate_mock_pdb
                    pdb_content = generate_mock_pdb(binder.sequence)
                else:
                    result = call_nvidia_protein_api(
                        binder.sequence,
                        selected_model["id"],
                        st.session_state.api_key,
                        selected_model_name
                    )
                    
                    if result["status"] == "success":
                        pdb_content = extract_pdb_from_response(result["data"])
                    else:
                        st.error(f"❌ Prediction failed: {result['message']}")
                        session.update_stage_status(WorkflowStage.BINDER_PREDICTION, StageStatus.FAILED)
                        return
                
                binder.pdb_content = pdb_content
                binder.structure_predicted = True
                binder.model_used = selected_model_name
                
                session.update_stage_status(WorkflowStage.BINDER_PREDICTION, StageStatus.COMPLETED)
                st.success("Binder structure prediction completed!")
                st.rerun()
                
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                session.update_stage_status(WorkflowStage.BINDER_PREDICTION, StageStatus.FAILED)
    
    # Show structure if available
    if binder.pdb_content:
        st.subheader("Binder Structure")
        try:
            # Binder structure from RFDiffusion doesn't have pLDDT - use solid color
            html_content = create_3d_visualization(binder.pdb_content, color_by_plddt=False)
            st.components.v1.html(html_content, height=500)
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
        
        st.download_button(
            "Download Binder PDB",
            data=binder.pdb_content,
            file_name=f"binder_{session.session_id}.pdb",
            mime="chemical/x-pdb"
        )
    
    # Navigation
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back", use_container_width=True):
            session.advance_to_stage(WorkflowStage.BINDER_DESIGN)
            st.rerun()
    
    with col2:
        if st.button("Next: Analyze Complex →", type="primary", use_container_width=True):
            if binder.pdb_content:
                session.advance_to_stage(WorkflowStage.COMPLEX_ANALYSIS)
                st.rerun()
            else:
                st.error("❌ Please predict binder structure first")


def render_complex_analysis_stage():
    """Stage 5: Complex Analysis"""
    st.header("5️⃣ Complex Analysis")
    st.markdown("Analyze the binding interface and interaction quality")
    
    session = st.session_state.workflow_session
    target = session.target
    binder = session.binder
    complex_data = session.complex
    
    # Docking method selection
    st.subheader("Docking Method")
    docking_method = st.radio(
        "Select Method",
        ["Simple Overlay", "DiffDock (Coming Soon)"],
        horizontal=True,
        help="Simple overlay places structures as-is. DiffDock predicts optimal binding pose."
    )
    
    complex_data.docking_method = docking_method.lower().split()[0]
    
    # Analysis parameters
    with st.expander("⚙️ Analysis Parameters"):
        interface_cutoff = st.slider(
            "Interface Distance Cutoff (Å)",
            min_value=3.0,
            max_value=8.0,
            value=5.0,
            step=0.5,
            help="Maximum distance for residues to be considered at the interface"
        )
    
    # Analyze button
    if st.button("Analyze Binding Interface", type="primary"):
        session.update_stage_status(WorkflowStage.COMPLEX_ANALYSIS, StageStatus.IN_PROGRESS)
        
        with st.spinner("Analyzing binding interface..."):
            try:
                # Parse structures
                target_atoms = parse_pdb_content(target.pdb_content)
                binder_atoms = parse_pdb_content(binder.pdb_content)
                
                # Find interface
                interface_data = find_interface_residues(
                    target_atoms, 
                    binder_atoms, 
                    cutoff=interface_cutoff
                )
                
                # Assess quality
                quality_data = assess_binding_quality(interface_data)
                
                # Combine structures
                combined_pdb = combine_pdbs(target.pdb_content, binder.pdb_content)
                
                # Store results
                complex_data.complex_pdb = combined_pdb
                complex_data.interface_residues_target = interface_data['interface_residues_target']
                complex_data.interface_residues_binder = interface_data['interface_residues_binder']
                complex_data.num_contacts = interface_data['num_contacts']
                complex_data.avg_distance = interface_data['avg_distance']
                complex_data.min_distance = interface_data['min_distance']
                complex_data.quality_score = quality_data['quality_score']
                complex_data.quality_grade = quality_data['grade']
                complex_data.feedback = quality_data['feedback']
                
                session.update_stage_status(WorkflowStage.COMPLEX_ANALYSIS, StageStatus.COMPLETED)
                st.success("✅ Analysis completed!")
                st.rerun()
                
            except Exception as e:
                st.error(f"❌ Analysis failed: {str(e)}")
                session.update_stage_status(WorkflowStage.COMPLEX_ANALYSIS, StageStatus.FAILED)
    
    # Show results if available
    if complex_data.complex_pdb:
        render_analysis_results(session)
    
    # Navigation
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back", use_container_width=True):
            session.advance_to_stage(WorkflowStage.BINDER_PREDICTION)
            st.rerun()
    
    with col2:
        if st.button("View Results →", type="primary", use_container_width=True):
            if complex_data.complex_pdb:
                session.advance_to_stage(WorkflowStage.RESULTS)
                st.rerun()
            else:
                st.error("❌ Please run analysis first")


def render_analysis_results(session: WorkflowSession):
    """Render analysis results"""
    complex_data = session.complex
    
    st.subheader("Analysis Results")
    
    # Quality metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("pLDDT Score", f"{complex_data.plddt_score or complex_data.quality_score or 0:.1f}")
    with col2:
        st.metric("Grade", complex_data.quality_grade)
    with col3:
        st.metric("Interface Contacts", complex_data.num_contacts)
    with col4:
        st.metric("Avg Distance", f"{complex_data.avg_distance:.2f} Å")
    
    # Feedback
    st.markdown("**Feedback:**")
    for feedback in complex_data.feedback:
        st.markdown(f"- {feedback}")
    
    # Complex visualization
    st.subheader("🧬 Complex Structure")
    try:
        html_content = create_3d_visualization(complex_data.complex_pdb)
        st.components.v1.html(html_content, height=600)
    except Exception as e:
        st.error(f"Visualization error: {str(e)}")
    
    # Download
    st.download_button(
        "Download Complex PDB",
        data=complex_data.complex_pdb,
        file_name=f"complex_{session.session_id}.pdb",
        mime="chemical/x-pdb"
    )


def render_results_stage():
    """Stage 6: Final Results"""
    st.header("6️⃣ Final Results")
    st.markdown("Complete analysis of your binding protein design")
    
    session = st.session_state.workflow_session
    
    # Summary dashboard
    render_results_dashboard(session)
    
    # Navigation
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back to Analysis", key="results_back_to_analysis", width="stretch"):
            session.advance_to_stage(WorkflowStage.COMPLEX_ANALYSIS)
            st.rerun()
    
    with col2:
        if st.button("🔄 Start New Design", key="results_start_new_design", width="stretch"):
            st.session_state.workflow_session = WorkflowSession.create_new()
            st.rerun()


def render_results_dashboard(session: WorkflowSession):
    """Render complete results dashboard - comprehensive summary of entire workflow"""
    target = session.target
    binder = session.binder
    complex_data = session.complex
    
    # Header with overall status
    quality_score = complex_data.quality_score or 0
    recommendation_color = (
        "[*]" if quality_score >= 70 else
        "🟡" if quality_score >= 50 else "🔴"
    )
    
    plddt_display = complex_data.plddt_score or quality_score
    st.markdown(f"""
    <div style="background: linear-gradient(135deg, #76B900 0%, #00D4AA 100%); padding: 30px; border-radius: 15px; color: white; margin-bottom: 30px;">
        <h1 style="margin: 0; font-size: 36px;">✅ Workflow Complete</h1>
        <p style="margin: 10px 0 0 0; font-size: 18px; opacity: 0.9;">
            Project: {session.project_name} | Status: {recommendation_color} pLDDT Score: {plddt_display:.1f}
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # === SECTION 1: Project Overview ===
    st.markdown("## 📋 Project Overview")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Session ID", session.session_id[:8])
        st.caption(f"Created: {session.created_at[:10]}")
    
    with col2:
        st.metric("Total Stages", "6")
        completed_count = sum(1 for status in session.stage_statuses.values() if status == StageStatus.COMPLETED.value)
        st.caption(f"Completed: {completed_count}/6")
    
    with col3:
        plddt_val = complex_data.plddt_score or quality_score
        st.metric("pLDDT Score", f"{plddt_val:.1f}")
        st.caption(f"Grade: {complex_data.quality_grade or 'N/A'}")
    
    st.markdown("---")
    
    # === SECTION 2: Target Protein Summary ===
    st.markdown("## 🎯 Target Protein Summary")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown(f"""
        **Input Details:**
        - **Input Type:** {target.input_type or 'N/A'}
        - **Sequence Length:** {len(target.sequence) if target.sequence else 0} amino acids
        - **Structure Prediction Model:** {target.model_used or 'Uploaded'}
        - **Structure Confidence (pLDDT):** {f"{target.confidence_avg:.2f}" if target.confidence_avg else 'N/A'}
        """)
        
        if target.sequence:
            with st.expander("View Target Sequence"):
                st.code(target.sequence, language=None)
    
    with col2:
        st.markdown("**Status:**")
        st.success("✅ Structure Ready")
        if target.structure_file_path:
            st.caption(f"File: {target.structure_file_path.split('/')[-1]}")
    
    st.markdown("---")
    
    # === SECTION 3: Binder Design Summary ===
    st.markdown("## 🧬 Binder Design Summary")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        binder_seq = binder.sequence or (binder.mpnn_sequences[binder.selected_sequence_idx] if binder.mpnn_sequences else None)
        
        st.markdown(f"""
        **Design Details:**
        - **Scaffold Generation:** RFDiffusion
        - **Sequence Design:** ProteinMPNN
        - **Selected Sequence:** #{binder.selected_sequence_idx + 1 if binder.mpnn_sequences else 'N/A'}
        - **Sequence Length:** {len(binder_seq) if binder_seq else 0} amino acids
        - **Total Sequences Generated:** {len(binder.mpnn_sequences) if binder.mpnn_sequences else 0}
        """)
        
        if binder.mpnn_scores and binder.selected_sequence_idx < len(binder.mpnn_scores):
            st.markdown(f"- **MPNN Score:** {binder.mpnn_scores[binder.selected_sequence_idx]:.3f}")
        
        if binder_seq:
            with st.expander("View Binder Sequence"):
                st.code(binder_seq, language=None)
    
    with col2:
        st.markdown("**Status:**")
        st.success("✅ Design Complete")
        if binder.scaffold_file_path:
            st.caption(f"Scaffold: {binder.scaffold_file_path.split('/')[-1]}")
    
    st.markdown("---")
    
    # === SECTION 4: Complex Analysis Summary ===
    st.markdown("## 🔗 Complex Binding Analysis")
    
    # Show only pLDDT score
    plddt_val = complex_data.plddt_score or quality_score
    st.metric("pLDDT Score", f"{plddt_val:.1f}")
    
    # Assessment feedback
    st.markdown("### 📊 Binding Assessment")
    feedback = complex_data.feedback or ["No binding analysis available yet"]
    for item in feedback:
        if "Strong" in item or "Excellent" in item or "Good" in item:
            st.success(f"✅ {item}")
        elif "Weak" in item or "Poor" in item or "Warning" in item:
            st.warning(f"⚠️ {item}")
        else:
            st.info(f"ℹ️ {item}")
    
    # Interface residues detail
    with st.expander("🔬 View Interface Residues Detail"):
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Target Interface Residues:**")
            target_interface = complex_data.interface_residues_target or []
            if target_interface:
                st.text(", ".join(map(str, target_interface)))
                st.caption(f"Total: {len(target_interface)} residues")
            else:
                st.text("N/A")
        
        with col2:
            st.markdown("**Binder Interface Residues:**")
            binder_interface = complex_data.interface_residues_binder or []
            if binder_interface:
                st.text(", ".join(map(str, binder_interface)))
                st.caption(f"Total: {len(binder_interface)} residues")
            else:
                st.text("N/A")
    
    st.markdown("---")
    
    # === SECTION 5: 3D Visualization ===
    st.markdown("## 🔬 3D Structure Visualization")
    
    if complex_data.complex_pdb:
        tab1, tab2, tab3 = st.tabs(["🔗 Complex", "🎯 Target", "🧬 Binder"])
        
        with tab1:
            st.markdown("**Target-Binder Complex**")
            try:
                from frontend.app_v2 import render_protein_viewer
                render_protein_viewer(complex_data.complex_pdb, key="results_complex_viewer")
            except:
                st.info("💡 3D visualization not available. Download the PDB file to view in external software.")
        
        with tab2:
            st.markdown("**Target Protein Structure**")
            if target.pdb_content:
                try:
                    from frontend.app_v2 import render_protein_viewer
                    render_protein_viewer(target.pdb_content, key="results_target_viewer")
                except:
                    st.info("💡 Download PDB to view externally.")
            else:
                st.warning("Target structure not available")
        
        with tab3:
            st.markdown("**Binder Protein Structure**")
            binder_pdb = binder.pdb_content or binder.scaffold_pdb
            if binder_pdb:
                try:
                    from frontend.app_v2 import render_protein_viewer
                    render_protein_viewer(binder_pdb, key="results_binder_viewer")
                except:
                    st.info("💡 Download PDB to view externally.")
            else:
                st.warning("Binder structure not available")
    else:
        st.info("Complex structure not available for visualization")
    
    st.markdown("---")
    
    # === SECTION 6: Export & Download ===
    st.markdown("## 📤 Export & Download Results")
    
    st.markdown("### Structure Files")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if target.pdb_content:
            st.download_button(
                "⬇️ Target PDB",
                data=target.pdb_content.encode('utf-8') if isinstance(target.pdb_content, str) else target.pdb_content,
                file_name=f"{session.project_name}_target.pdb",
                mime="chemical/x-pdb",
                use_container_width=True
            )
        else:
            st.button("⬇️ Target PDB", disabled=True, use_container_width=True)
    
    with col2:
        binder_pdb = binder.pdb_content or binder.scaffold_pdb
        if binder_pdb:
            st.download_button(
                "⬇️ Binder PDB",
                data=binder_pdb.encode('utf-8') if isinstance(binder_pdb, str) else binder_pdb,
                file_name=f"{session.project_name}_binder.pdb",
                mime="chemical/x-pdb",
                use_container_width=True
            )
        else:
            st.button("⬇️ Binder PDB", disabled=True, use_container_width=True)
    
    with col3:
        if complex_data.complex_pdb:
            st.download_button(
                "⬇️ Complex PDB",
                data=complex_data.complex_pdb.encode('utf-8') if isinstance(complex_data.complex_pdb, str) else complex_data.complex_pdb,
                file_name=f"{session.project_name}_complex.pdb",
                mime="chemical/x-pdb",
                use_container_width=True
            )
        else:
            st.button("⬇️ Complex PDB", disabled=True, use_container_width=True)
    
    st.markdown("### Sequence Files")
    col1, col2 = st.columns(2)
    
    with col1:
        if binder.mpnn_fasta_content:
            st.download_button(
                "⬇️ All MPNN Sequences (FASTA)",
                data=binder.mpnn_fasta_content,
                file_name=f"{session.project_name}_sequences.fa",
                mime="text/plain",
                use_container_width=True
            )
        else:
            st.button("⬇️ All MPNN Sequences", disabled=True, use_container_width=True)
    
    with col2:
        export_data = session.to_json()
        st.download_button(
            "⬇️ Complete Project (JSON)",
            data=export_data,
            file_name=f"{session.project_name}_complete.json",
            mime="application/json",
            use_container_width=True,
            type="primary"
        )
    
    st.markdown("---")
    
    # === SECTION 7: Next Steps ===
    st.markdown("## 🚀 Next Steps")
    
    # Use a more balanced 3-column layout
    col_rec, col_spacer, col_actions = st.columns([1, 0.1, 1])
    
    with col_rec:
        st.markdown("""
        **💡 Recommendations:**
        - Validate binding experimentally (e.g., SPR, ITC)
        - Optimize binder sequence for expression
        - Test different MPNN sequences
        - Refine complex with molecular dynamics
        """)
    
    with col_actions:
        st.markdown("**🔧 Actions:**")
        
        # Stack buttons vertically with consistent spacing
        if st.button("🔄 Design New Binder", use_container_width=True):
            session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
            st.rerun()
        
        if st.button("🆕 Start New Project", use_container_width=True):
            st.session_state.workflow_session = WorkflowSession.create_new()
            st.rerun()
    
    # Add navigation buttons at the bottom in a clean row
    st.markdown("---")
    nav_col1, nav_col2 = st.columns(2)
    
    with nav_col1:
        if st.button("← Back to Analysis", use_container_width=True):
            session.advance_to_stage(WorkflowStage.COMPLEX_PREDICTION)
            st.rerun()
    
    with nav_col2:
        if st.button("🔄 Start New Design", use_container_width=True, type="primary"):
            session.advance_to_stage(WorkflowStage.TARGET_INPUT)
            st.rerun()


def render_examples_gallery():
    """Render examples gallery page"""
    # Header
    st.markdown("""
    <div style="text-align: center; padding: 20px 0;">
        <h1>🧬 PROTEIN BINDER EXAMPLES GALLERY</h1>
        <p style="font-size: 18px; color: #6c757d; margin-top: 10px;">
            Explore computationally designed protein binders targeting key therapeutic targets
        </p>
    </div>
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
    
    # Info banner
    st.info("💡 **Tip:** Click 'View 3D Structure' on any example to interact - rotate, zoom, and inspect residues!")
    
    # Load examples directory
    examples_dir = Path(__file__).parent.parent / "examples_examples"
    
    if not examples_dir.exists():
        st.error(f"Examples directory not found: {examples_dir}")
        if st.button("🔙 Back to Workflow", type="primary"):
            st.session_state.show_examples = False
            st.rerun()
        return
    
    # Render examples
    for filename, metadata in EXAMPLES.items():
        # Example card
        st.markdown(f"""
        <div style="background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); 
                    border-left: 4px solid #76B900; border-radius: 8px; 
                    padding: 20px; margin-bottom: 20px; 
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
            <div style="color: #1A1A1A; font-size: 24px; font-weight: 700; margin-bottom: 10px;">
                {metadata['icon']} {metadata['name']}
            </div>
            <div style="color: #495057; font-size: 14px; line-height: 1.6; margin-bottom: 15px;">
                <strong>Target:</strong> {metadata['target']}<br>
                {metadata['description']}
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Metrics row
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Use Case", metadata['use_case'].split()[0])
        with col2:
            st.metric("Confidence", metadata['confidence'])
        with col3:
            st.metric("Applications", len(metadata['applications']))
        
        file_path = examples_dir / filename
        if file_path.exists():
            # Load binary PDB - handle both ZIP and plain text
            try:
                with open(file_path, 'rb') as f:
                    file_bytes = f.read()
                
                # Check if it's a ZIP file
                import zipfile
                import io
                
                pdb_content = None
                try:
                    # Try to read as ZIP first (EGFR and Insulin examples)
                    with zipfile.ZipFile(io.BytesIO(file_bytes)) as zf:
                        # Extract the first file which contains PDB text
                        first_file = zf.namelist()[0]
                        pdb_content = zf.read(first_file).decode('utf-8', errors='ignore')
                
                except zipfile.BadZipFile:
                    # Not a ZIP, try as plain text (COVID_SPIKE, KRAS_G12D, pdl_1)
                    pdb_content = file_bytes.decode('utf-8', errors='ignore')
                
                if not pdb_content or 'ATOM' not in pdb_content:
                    st.error(f"Could not extract valid PDB data from {filename}")
                    with col4:
                        st.metric("Atoms", "N/A")
                    continue
                
                with col4:
                    atom_count = pdb_content.count('ATOM')
                    st.metric("Atoms", f"{atom_count:,}")
                
                # 3D Visualization
                with st.expander("🔬 View 3D Structure", expanded=False):
                    try:
                        html_content = create_3d_visualization(pdb_content, color_by_plddt=True)
                        components.html(html_content, height=600, scrolling=False)
                    except Exception as e:
                        st.error(f"Visualization error: {str(e)}")
                
                # Applications
                with st.expander("📋 Potential Applications"):
                    for app in metadata['applications']:
                        st.markdown(f"• **{app}**")
                        
            except Exception as e:
                st.error(f"Error loading {filename}: {str(e)}")
        else:
            st.warning(f"File not found: {filename}")
        
        st.markdown("<hr style='margin: 30px 0; border: 1px solid #e9ecef;'>", unsafe_allow_html=True)
    
    # Back button
    if st.button("🔙 Back to Workflow", type="primary", use_container_width=True):
        st.session_state.show_examples = False
        st.rerun()


def main():
    """Main application"""
    initialize_session_state()
    
    # Check for navigation from pipeline via query params
    query_params = st.query_params
    if 'stage' in query_params:
        stage_map = {
            'target_input': WorkflowStage.TARGET_INPUT,
            'target_prediction': WorkflowStage.TARGET_PREDICTION,
            'binder_scaffold': WorkflowStage.BINDER_SCAFFOLD_DESIGN,
            'binder_sequence': WorkflowStage.BINDER_SEQUENCE_DESIGN,
            'complex_prediction': WorkflowStage.COMPLEX_PREDICTION,
            'results': WorkflowStage.RESULTS
        }
        requested_stage = query_params['stage']
        if requested_stage in stage_map:
            st.session_state.workflow_session.advance_to_stage(stage_map[requested_stage])
            st.session_state.show_pipeline = False
            # Clear the query param
            st.query_params.clear()
            st.rerun()
    
    # Sidebar
    render_sidebar()
    
    # Show examples gallery if button was clicked AND user hasn't navigated to a different stage
    # If user clicks a stage nav button while in gallery, it should close the gallery and show that stage
    session = st.session_state.workflow_session
    if st.session_state.get('show_examples', False) and not st.session_state.get('_stage_nav_clicked', False):
        render_examples_gallery()
        return
    elif st.session_state.get('_stage_nav_clicked', False):
        # User clicked a stage navigation button, close gallery
        st.session_state.show_examples = False
        st.session_state._stage_nav_clicked = False
    
    # Show pipeline if button was clicked
    if st.session_state.get('show_pipeline', False):
        # Read and embed the HTML file directly using st.components.v1.html()
        pipeline_path = Path(__file__).parent / "pipeline.html"
        if pipeline_path.exists():
            with open(pipeline_path, 'r', encoding='utf-8') as f:
                pipeline_html = f.read()
            st.components.v1.html(pipeline_html, height=2500, scrolling=True)
        else:
            st.error("Pipeline diagram not found.")
        
        # Back button
        if st.button("🔙 Back to Workflow", use_container_width=True, type="primary"):
            st.session_state.show_pipeline = False
            st.rerun()
        return
    
    # NVIDIA Branded Header
    st.markdown("""
    <div class="nvidia-header">
        <h1>🧬 Protein Binding Design Workflow</h1>
        <p>AlphaFold2 → RFDiffusion → ProteinMPNN → AlphaFold Multimer</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Progress stepper
    render_progress_stepper()
    
    st.markdown("---")
    
    # Render current stage
    session = st.session_state.workflow_session
    current_stage = session.current_stage
    
    if current_stage == WorkflowStage.TARGET_INPUT:
        render_target_input_stage()
    elif current_stage == WorkflowStage.TARGET_PREDICTION:
        render_target_prediction_stage()
    elif current_stage == WorkflowStage.BINDER_SCAFFOLD_DESIGN:
        render_binder_scaffold_stage()
    elif current_stage == WorkflowStage.BINDER_SEQUENCE_DESIGN:
        render_binder_sequence_stage()
    elif current_stage == WorkflowStage.COMPLEX_PREDICTION:
        render_complex_prediction_stage()
    elif current_stage == WorkflowStage.RESULTS:
        render_results_stage()
    
    # NVIDIA Footer
    st.markdown("---")
    
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col2:
        st.markdown("""
        <div style="text-align: center; color: #666666; padding: 20px;">
            <p style="font-size: 14px; margin-bottom: 10px;">
                <strong>Powered by NVIDIA AI</strong><br>
                This workflow uses NVIDIA's state-of-the-art protein folding models for structure prediction.
            </p>
            <p style="font-size: 12px; color: #999999;">
                Results are computational predictions and should be validated experimentally.<br>
                © 2025 NVIDIA Corporation. All rights reserved.
            </p>
        </div>
        """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
