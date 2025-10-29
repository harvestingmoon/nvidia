"""
Binding Protein Design Workflow Application
Multi-step workflow for designing and analyzing protein binders
ESMFold → RFDiffusion → ProteinMPNN → DiffDock
"""

import streamlit as st
import time
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
import plotly.graph_objects as go
import plotly.express as px

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from workflow.workflow_state import (
    WorkflowSession, WorkflowStage, StageStatus, 
    WorkflowValidator, TargetProteinData, BinderProteinData
)
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
    }
    
    .stage-active {
        color: #00D4AA;
        font-weight: 700;
        font-size: 16px;
        animation: pulse 2s infinite;
    }
    
    .stage-pending {
        color: #666666;
        font-weight: 500;
        font-size: 16px;
    }
    
    .stage-failed {
        color: #FF3838;
        font-weight: 700;
        font-size: 16px;
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
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);
        color: white;
        border: none;
        border-radius: 6px;
        padding: 10px 24px;
        font-weight: 600;
        font-size: 16px;
        transition: all 0.3s;
        box-shadow: 0 2px 8px rgba(118, 185, 0, 0.3);
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #5A8F00 0%, #76B900 100%);
        box-shadow: 0 4px 12px rgba(118, 185, 0, 0.4);
        transform: translateY(-1px);
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


def initialize_session_state():
    """Initialize session state variables"""
    if 'workflow_session' not in st.session_state:
        st.session_state.workflow_session = WorkflowSession.create_new()
    
    if 'api_key' not in st.session_state:
        st.session_state.api_key = None
    
    if 'demo_mode' not in st.session_state:
        st.session_state.demo_mode = False


def render_progress_stepper():
    """Render workflow progress stepper"""
    session = st.session_state.workflow_session
    
    stages = [
        ("1️⃣ Target", WorkflowStage.TARGET_INPUT),
        ("2️⃣ Predict", WorkflowStage.TARGET_PREDICTION),
        ("3️⃣ Binder", WorkflowStage.BINDER_DESIGN),
        ("4️⃣ Predict", WorkflowStage.BINDER_PREDICTION),
        ("5️⃣ Analyze", WorkflowStage.COMPLEX_ANALYSIS),
        ("6️⃣ Results", WorkflowStage.RESULTS)
    ]
    
    st.markdown("### 🔬 Workflow Progress")
    
    cols = st.columns(len(stages))
    
    for idx, (col, (label, stage)) in enumerate(zip(cols, stages)):
        with col:
            status = session.stage_statuses.get(stage.value, "not_started")
            
            if status == StageStatus.COMPLETED.value:
                st.markdown(f'<p class="stage-complete">✅ {label}</p>', unsafe_allow_html=True)
            elif status == StageStatus.IN_PROGRESS.value:
                st.markdown(f'<p class="stage-active">🔄 {label}</p>', unsafe_allow_html=True)
            elif status == StageStatus.FAILED.value:
                st.markdown(f'<p style="color: #dc3545; font-weight: bold;">❌ {label}</p>', unsafe_allow_html=True)
            else:
                st.markdown(f'<p class="stage-pending">⭕ {label}</p>', unsafe_allow_html=True)
    
    st.progress(calculate_overall_progress(session))


def calculate_overall_progress(session: WorkflowSession) -> float:
    """Calculate overall workflow progress"""
    total_stages = len(WorkflowStage)
    completed = sum(1 for status in session.stage_statuses.values() 
                   if status == StageStatus.COMPLETED.value)
    return completed / total_stages


def render_sidebar():
    """Render sidebar with session management"""
    # NVIDIA Logo in Sidebar
    try:
        st.sidebar.image("image/nvidia.jpg", use_container_width=True)
    except:
        st.sidebar.markdown("# NVIDIA")
    
    st.sidebar.title("🧬 Binding Workflow")
    st.sidebar.markdown("---")
    
    session = st.session_state.workflow_session
    
    # Project info
    st.sidebar.subheader("📁 Project")
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
    
    st.sidebar.markdown("---")
    
    # API Configuration
    st.sidebar.subheader("⚙️ Configuration")
    
    demo_mode = st.sidebar.checkbox(
        "Demo Mode (No API)",
        value=st.session_state.demo_mode,
        help="Use demo mode to test the interface without API calls"
    )
    st.session_state.demo_mode = demo_mode
    
    if not demo_mode:
        api_key = st.sidebar.text_input(
            "NVIDIA API Key",
            type="password",
            value=st.session_state.api_key or "",
            help="Get your API key from NVIDIA NGC"
        )
        if api_key:
            st.session_state.api_key = api_key
    
    st.sidebar.markdown("---")
    
    # Session management
    st.sidebar.subheader("💾 Session")
    
    col1, col2 = st.sidebar.columns(2)
    
    with col1:
        if st.button("💾 Save", use_container_width=True):
            save_session()
    
    with col2:
        if st.button("📂 Load", use_container_width=True):
            load_session()
    
    if st.sidebar.button("🔄 Reset Workflow", use_container_width=True):
        if st.sidebar.checkbox("Confirm reset?"):
            st.session_state.workflow_session = WorkflowSession.create_new()
            st.rerun()
    
    # Export options
    with st.sidebar.expander("📤 Export Options"):
        if st.button("Export as JSON"):
            export_json = session.to_json()
            st.download_button(
                "Download JSON",
                data=export_json,
                file_name=f"{session.project_name}_{session.session_id}.json",
                mime="application/json"
            )
    
    # Notes
    st.sidebar.markdown("---")
    st.sidebar.subheader("📝 Notes")
    notes = st.sidebar.text_area(
        "Project Notes",
        value=session.notes,
        height=100,
        key="notes_input"
    )
    if notes != session.notes:
        session.notes = notes


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


def render_target_input_stage():
    """Stage 1: Target Protein Input"""
    st.header("1️⃣ Target Protein Input")
    st.markdown("Define the protein you want to design a binder for")
    
    session = st.session_state.workflow_session
    target = session.target
    
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
    with st.expander("⚙️ Advanced: Specify Binding Site (Optional)"):
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
    
    session = st.session_state.workflow_session
    target = session.target
    
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
        # Model selection
        st.subheader("Structure Prediction")
        
        model_options = {name: model for name, model in PROTEIN_MODELS.items()}
        selected_model_name = st.selectbox(
            "Select Prediction Model",
            options=list(model_options.keys()),
            index=1,  # Default to OpenFold2
            help="OpenFold2 is faster, AlphaFold2 is more accurate"
        )
        
        selected_model = model_options[selected_model_name]
        
        # Predict button
        if st.button("🔬 Predict Structure", type="primary"):
            if not st.session_state.api_key and not st.session_state.demo_mode:
                st.error("❌ Please provide an API key or enable Demo Mode")
                return
            
            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.IN_PROGRESS)
            
            with st.spinner(f"Predicting structure with {selected_model_name}..."):
                try:
                    if st.session_state.demo_mode:
                        # Demo mode
                        time.sleep(2)
                        from frontend.app_v2 import generate_mock_pdb
                        pdb_content = generate_mock_pdb(target.sequence)
                    else:
                        # Real prediction
                        result = call_nvidia_protein_api(
                            target.sequence,
                            selected_model["id"],
                            st.session_state.api_key,
                            selected_model_name
                        )
                        
                        if result["status"] == "success":
                            pdb_content = extract_pdb_from_response(result["data"])
                        else:
                            st.error(f"❌ Prediction failed: {result['message']}")
                            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
                            return
                    
                    # Store results
                    target.pdb_content = pdb_content
                    target.structure_predicted = True
                    target.model_used = selected_model_name
                    
                    session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
                    st.success("🎉 Structure prediction completed!")
                    st.rerun()
                    
                except Exception as e:
                    st.error(f"❌ Error: {str(e)}")
                    session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
    
    # Show structure if available
    if target.pdb_content:
        st.subheader("Target Structure")
        try:
            html_content = create_3d_visualization(target.pdb_content)
            st.components.v1.html(html_content, height=500)
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
        
        # Download button
        st.download_button(
            "📥 Download Target PDB",
            data=target.pdb_content,
            file_name=f"target_{session.session_id}.pdb",
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
                session.advance_to_stage(WorkflowStage.BINDER_DESIGN)
                st.rerun()
            else:
                st.error("❌ Please predict target structure first")


def render_binder_design_stage():
    """Stage 3: Binder Design"""
    st.header("3️⃣ Binder Design")
    st.markdown("Design or input the binding protein sequence")
    
    session = st.session_state.workflow_session
    binder = session.binder
    
    # Design method selection
    design_method = st.radio(
        "Design Method",
        ["Manual Sequence", "RFDiffusion (Coming Soon)", "Template-based (Coming Soon)"],
        horizontal=True
    )
    
    binder.design_method = design_method.lower().split()[0]
    
    if design_method == "Manual Sequence":
        st.subheader("Enter Binder Sequence")
        
        # Sequence length guidance
        col1, col2 = st.columns([2, 1])
        
        with col1:
            sequence = st.text_area(
                "Binder Sequence",
                value=binder.sequence or "",
                placeholder="Enter designed binder sequence (50-150 residues recommended)",
                height=150
            )
        
        with col2:
            st.markdown("**💡 Design Tips:**")
            st.markdown("- Typical binder: 50-150 AA")
            st.markdown("- Include diverse residues")
            st.markdown("- Consider hydrophobic core")
            st.markdown("- Plan interface residues")
        
        if sequence:
            is_valid, result = validate_protein_sequence(sequence)
            if is_valid:
                binder.sequence = result
                st.success(f"✅ Valid binder sequence: {len(result)} amino acids")
                
                # Sequence analysis
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Length", len(result))
                with col2:
                    hydrophobic = sum(1 for aa in result if aa in 'AILMFVPW')
                    st.metric("Hydrophobic %", f"{hydrophobic/len(result)*100:.1f}%")
                with col3:
                    charged = sum(1 for aa in result if aa in 'DEKR')
                    st.metric("Charged %", f"{charged/len(result)*100:.1f}%")
            else:
                st.error(f"❌ {result}")
    
    else:
        st.info("🚧 This feature is coming soon! Use Manual Sequence for now.")
    
    # Navigation
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("← Back", use_container_width=True):
            session.advance_to_stage(WorkflowStage.TARGET_PREDICTION)
            st.rerun()
    
    with col2:
        if st.button("Next: Predict Binder →", type="primary", use_container_width=True):
            if binder.sequence:
                session.advance_to_stage(WorkflowStage.BINDER_PREDICTION)
                session.update_stage_status(WorkflowStage.BINDER_DESIGN, StageStatus.COMPLETED)
                st.rerun()
            else:
                st.error("❌ Please enter a binder sequence")


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
    if st.button("🔬 Predict Binder Structure", type="primary"):
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
                st.success("🎉 Binder structure prediction completed!")
                st.rerun()
                
            except Exception as e:
                st.error(f"❌ Error: {str(e)}")
                session.update_stage_status(WorkflowStage.BINDER_PREDICTION, StageStatus.FAILED)
    
    # Show structure if available
    if binder.pdb_content:
        st.subheader("Binder Structure")
        try:
            html_content = create_3d_visualization(binder.pdb_content)
            st.components.v1.html(html_content, height=500)
        except Exception as e:
            st.error(f"Visualization error: {str(e)}")
        
        st.download_button(
            "📥 Download Binder PDB",
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
    if st.button("🔬 Analyze Binding Interface", type="primary"):
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
    
    st.subheader("📊 Analysis Results")
    
    # Quality metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Quality Score", f"{complex_data.quality_score}/100")
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
        "📥 Download Complex PDB",
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
        if st.button("← Back to Analysis", use_container_width=True):
            session.advance_to_stage(WorkflowStage.COMPLEX_ANALYSIS)
            st.rerun()
    
    with col2:
        if st.button("🔄 Start New Design", use_container_width=True):
            st.session_state.workflow_session = WorkflowSession.create_new()
            st.rerun()


def render_results_dashboard(session: WorkflowSession):
    """Render complete results dashboard"""
    target = session.target
    binder = session.binder
    complex_data = session.complex
    
    # Overall summary
    st.subheader("📋 Project Summary")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Target Protein:**")
        st.markdown(f"- Length: {len(target.sequence)} AA")
        st.markdown(f"- Model: {target.model_used or 'Uploaded'}")
        st.markdown(f"- Input: {target.input_type}")
    
    with col2:
        st.markdown("**Binder Protein:**")
        st.markdown(f"- Length: {len(binder.sequence)} AA")
        st.markdown(f"- Model: {binder.model_used}")
        st.markdown(f"- Design: {binder.design_method}")
    
    st.markdown("---")
    
    # Quality assessment
    st.subheader("🎯 Binding Quality Assessment")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            "Overall Quality Score",
            f"{complex_data.quality_score}/100",
            delta=None
        )
    
    with col2:
        st.metric("Grade", complex_data.quality_grade)
    
    with col3:
        recommendation_color = (
            "🟢" if complex_data.quality_score >= 70 else
            "🟡" if complex_data.quality_score >= 50 else "🔴"
        )
        st.markdown(f"**Status:** {recommendation_color}")
    
    # Detailed feedback
    st.markdown("**Assessment:**")
    for item in complex_data.feedback:
        st.markdown(f"- {item}")
    
    st.markdown("---")
    
    # Interface details
    st.subheader("🔗 Interface Details")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Interface Contacts", complex_data.num_contacts)
    with col2:
        st.metric("Average Distance", f"{complex_data.avg_distance:.2f} Å")
    with col3:
        st.metric("Minimum Distance", f"{complex_data.min_distance:.2f} Å")
    
    # Interface residues
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Target Interface Residues:**")
        st.text(", ".join(map(str, complex_data.interface_residues_target[:20])))
        if len(complex_data.interface_residues_target) > 20:
            st.caption(f"... and {len(complex_data.interface_residues_target) - 20} more")
    
    with col2:
        st.markdown("**Binder Interface Residues:**")
        st.text(", ".join(map(str, complex_data.interface_residues_binder[:20])))
        if len(complex_data.interface_residues_binder) > 20:
            st.caption(f"... and {len(complex_data.interface_residues_binder) - 20} more")
    
    st.markdown("---")
    
    # Export all results
    st.subheader("📤 Export Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.download_button(
            "Download Target PDB",
            data=target.pdb_content,
            file_name=f"target_{session.session_id}.pdb",
            mime="chemical/x-pdb"
        )
    
    with col2:
        st.download_button(
            "Download Binder PDB",
            data=binder.pdb_content,
            file_name=f"binder_{session.session_id}.pdb",
            mime="chemical/x-pdb"
        )
    
    with col3:
        st.download_button(
            "Download Complex PDB",
            data=complex_data.complex_pdb,
            file_name=f"complex_{session.session_id}.pdb",
            mime="chemical/x-pdb"
        )
    
    # Full session export
    if st.button("📦 Export Complete Project"):
        export_data = session.to_json()
        st.download_button(
            "Download Project JSON",
            data=export_data,
            file_name=f"{session.project_name}_{session.session_id}.json",
            mime="application/json"
        )


def main():
    """Main application"""
    initialize_session_state()
    
    # Sidebar
    render_sidebar()
    
    # NVIDIA Branded Header
    st.markdown("""
    <div class="nvidia-header">
        <h1>🧬 NVIDIA Protein Binding Design Workflow</h1>
        <p>ESMFold → RFDiffusion → ProteinMPNN → DiffDock</p>
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
    elif current_stage == WorkflowStage.BINDER_DESIGN:
        render_binder_design_stage()
    elif current_stage == WorkflowStage.BINDER_PREDICTION:
        render_binder_prediction_stage()
    elif current_stage == WorkflowStage.COMPLEX_ANALYSIS:
        render_complex_analysis_stage()
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
