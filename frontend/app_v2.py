"""
Protein Structure Prediction with NVIDIA Cloud Functions
A Streamlit application for predicting protein 3D structures using NVIDIA's protein folding models
"""

import streamlit as st
import requests
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import py3Dmol
import streamlit.components.v1 as components
import re
import time
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.protein_models import PROTEIN_MODELS

# Configure Streamlit page
st.set_page_config(
    page_title="NVIDIA Protein Structure Prediction",
    page_icon="🧬",
    layout="wide"
)

# NVIDIA Theme CSS
st.markdown("""
<style>
    /* NVIDIA Brand Colors */
    :root {
        --nvidia-green: #76B900;
        --nvidia-dark: #1A1A1A;
    }
    
    /* Headers with NVIDIA styling */
    h1, h2, h3 {
        color: #1A1A1A !important;
        font-family: 'NVIDIA Sans', Arial, sans-serif;
    }
    
    h1 {
        border-bottom: 3px solid #76B900;
        padding-bottom: 10px;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #76B900 0%, #5A8F00 100%);
        color: white;
        border: none;
        border-radius: 6px;
        padding: 10px 24px;
        font-weight: 600;
        transition: all 0.3s;
    }
    
    .stButton > button:hover {
        background: linear-gradient(135deg, #5A8F00 0%, #76B900 100%);
        transform: translateY(-1px);
        box-shadow: 0 4px 12px rgba(118, 185, 0, 0.4);
    }
    
    /* Success messages */
    .stSuccess {
        background-color: rgba(118, 185, 0, 0.1);
        border-left: 4px solid #76B900;
    }
    
    /* Metrics */
    [data-testid="stMetricValue"] {
        color: #76B900;
        font-weight: 700;
    }
    
    /* Download button */
    .stDownloadButton > button {
        background: linear-gradient(135deg, #00D4AA 0%, #00A88A 100%);
        color: white;
        font-weight: 600;
    }
    
    /* Sidebar */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1A1A1A 0%, #2D2D2D 100%);
    }
    
    [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2 {
        color: #76B900 !important;
    }
    
    [data-testid="stSidebar"] p, [data-testid="stSidebar"] label {
        color: #E5E5E5 !important;
    }
</style>
""", unsafe_allow_html=True)

def validate_protein_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate if the input is a valid amino acid sequence
    
    Args:
        sequence (str): Input protein sequence
        
    Returns:
        tuple: (is_valid, cleaned_sequence_or_error_message)
    """
    if not sequence or len(sequence.strip()) == 0:
        return False, "Please enter a protein sequence"
    
    # Remove whitespace and convert to uppercase
    clean_sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Check if sequence contains only valid amino acid codes
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(clean_sequence) - valid_amino_acids
    
    if invalid_chars:
        return False, f"Invalid amino acid characters found: {', '.join(invalid_chars)}"
    
    # Check minimum length
    if len(clean_sequence) < 10:
        return False, "Sequence too short. Please enter at least 10 amino acids"
    
    # Check maximum length
    if len(clean_sequence) > 2000:
        return False, "Sequence too long. Please enter less than 2000 amino acids"
    
    return True, clean_sequence

def call_nvidia_protein_api(sequence: str, model_id: str, api_key: str, model_name: str = "Unknown") -> Dict[str, Any]:
    """
    Call NVIDIA Cloud Functions API for protein structure prediction
    
    Args:
        sequence (str): Protein amino acid sequence
        model_id (str): Model function ID
        api_key (str): NVIDIA API key
        model_name (str): Name of the model (default is "Unknown")
        
    Returns:
        dict: API response with status and data
    """
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "NVCF-POLL-SECONDS": "300"
    }
    
    # Different payload formats to try
    payloads_to_try = []
    
    # AlphaFold2 specific format (NVIDIA Health API) - ONLY this format for AlphaFold2
    if "alphafold" in model_name.lower():
        alphafold2_payload = {
            "sequence": sequence,
            "algorithm": "mmseqs2",
            "e_value": 0.0001,
            "iterations": 1,
            "databases": ["small_bfd"],
            "relax_prediction": False,
            "skip_template_search": True
        }
        payloads_to_try.append(alphafold2_payload)
    else:
        # Generic formats for other models
        payloads_to_try.extend([
            {"sequence": sequence},
            {"input": sequence},
            {"protein_sequence": sequence},
            {"sequences": [sequence]},
            {"data": {"sequence": sequence}},
            {
                "sequence": sequence,
                "num_recycles": 3,
                "max_templates": 4
            }
        ])
    
    # Use different endpoints based on model
    if "alphafold" in model_name.lower():
        endpoint = "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2"
        status_endpoint = "https://health.api.nvidia.com/v1/status"
    else:
        endpoint = f"https://api.nvcf.nvidia.com/v2/nvcf/pexec/functions/{model_id}"
        status_endpoint = None
    
    # Set longer timeout for AlphaFold2 models
    timeout_seconds = 600 if "alphafold" in model_name.lower() else 300
    
    for i, payload in enumerate(payloads_to_try):
        try:
            st.info(f"Trying payload format {i+1}...")
            
            # Show specific info for AlphaFold2
            if "alphafold" in model_name.lower():
                st.warning("⏳ AlphaFold2 may take 5-10 minutes for structure prediction. Please be patient...")
                if i == 0:
                    st.info("🧬 Using NVIDIA Health API format with MSA search...")
            
            response = requests.post(
                endpoint,
                headers=headers,
                json=payload,
                timeout=timeout_seconds
            )
            
            if response.status_code == 200:
                return {"status": "success", "data": response.json()}
            elif response.status_code == 202:
                # Asynchronous processing
                if "alphafold" in model_name.lower() and status_endpoint:
                    # AlphaFold2 specific polling
                    req_id = response.headers.get("nvcf-reqid")
                    if req_id:
                        st.success(f"✅ AlphaFold2 request accepted! Request ID: {req_id}")
                        return poll_alphafold2_result(req_id, api_key, status_endpoint, model_name)
                else:
                    # Generic polling
                    result = response.json()
                    if "reqId" in result:
                        st.success(f"✅ Request accepted! Processing with {model_name}...")
                        return poll_for_result(result["reqId"], api_key, model_name)
                    else:
                        return {"status": "error", "message": f"Async request submitted but no ID received: {response.text}"}
            elif response.status_code == 504:
                st.warning(f"⏱️ Payload format {i+1} timed out (504). This is common with {model_name}. Trying next format...")
                continue
            else:
                error_msg = f"Status {response.status_code}"
                try:
                    error_detail = response.json().get("detail", response.text)
                    error_msg += f": {error_detail}"
                except:
                    error_msg += f": {response.text}"
                
                st.warning(f"❌ Payload format {i+1} failed with {error_msg}")
                continue
                
        except requests.exceptions.Timeout:
            st.warning(f"⏱️ Payload format {i+1} timed out after {timeout_seconds} seconds. Trying next format...")
            continue
        except requests.exceptions.RequestException as e:
            st.warning(f"❌ Payload format {i+1} failed with exception: {str(e)}")
            continue
    
    return {
        "status": "error", 
        "message": f"All payload formats failed for {model_name}. This could be due to:\n"
                   f"• Server overload (common with AlphaFold2)\n"
                   f"• Sequence too long or complex\n"
                   f"• API rate limiting\n"
                   f"• Model temporarily unavailable\n\n"
                   f"💡 Try again in a few minutes or switch to OpenFold2 which is typically faster."
    }

def poll_for_result(request_id: str, api_key: str, model_name: str = "Unknown", max_attempts: int = 120) -> Dict[str, Any]:
    """
    Poll for the result of an asynchronous request
    Extended timeout for AlphaFold2 models
    """
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    # Increase max attempts for AlphaFold2
    if "alphafold" in model_name.lower():
        max_attempts = 180  # 30 minutes at 10-second intervals
    
    progress_bar = st.progress(0)
    status_placeholder = st.empty()
    
    for attempt in range(max_attempts):
        try:
            progress = attempt / max_attempts
            progress_bar.progress(progress)
            
            minutes_elapsed = (attempt * 10) // 60
            seconds_elapsed = (attempt * 10) % 60
            
            status_placeholder.info(
                f"🔬 {model_name} is processing your protein structure...\n"
                f"Time elapsed: {minutes_elapsed}m {seconds_elapsed}s | "
                f"Attempt {attempt + 1}/{max_attempts}"
            )
            
            # Poll the status endpoint
            poll_response = requests.get(
                f"https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{request_id}",
                headers=headers,
                timeout=30
            )
            
            if poll_response.status_code == 200:
                result = poll_response.json()
                status = result.get("status", "").upper()
                
                if status == "COMPLETED":
                    # Get the final result
                    result_response = requests.get(
                        f"https://api.nvcf.nvidia.com/v2/nvcf/pexec/response/{request_id}",
                        headers=headers,
                        timeout=60
                    )
                    
                    if result_response.status_code == 200:
                        progress_bar.progress(1.0)
                        status_placeholder.success(f"🎉 {model_name} structure prediction completed!")
                        return {"status": "success", "data": result_response.json()}
                    else:
                        return {"status": "error", "message": f"Failed to get result: {result_response.status_code} - {result_response.text}"}
                
                elif status == "FAILED":
                    error_msg = result.get("error", "Unknown error")
                    return {"status": "error", "message": f"{model_name} prediction failed: {error_msg}"}
                
                elif status in ["PENDING", "IN_PROGRESS", "QUEUED"]:
                    if status == "IN_PROGRESS":
                        status_placeholder.info(f"🧬 {model_name} is actively processing your sequence... (Step {attempt + 1})")
                    elif status == "QUEUED":
                        status_placeholder.info(f"⏳ Your {model_name} request is queued... (Position in queue: Step {attempt + 1})")
                    
                    time.sleep(10)
                    continue
                
                else:
                    return {"status": "error", "message": f"Unknown status from {model_name}: {status}"}
            
            else:
                time.sleep(10)
                continue
                
        except requests.exceptions.RequestException as e:
            if attempt == max_attempts - 1:
                return {"status": "error", "message": f"{model_name} polling failed after {max_attempts} attempts: {str(e)}"}
            time.sleep(10)
            continue
    
    timeout_minutes = (max_attempts * 10) // 60
    return {
        "status": "error", 
        "message": f"⏱️ {model_name} prediction timed out after {timeout_minutes} minutes.\n\n"
                   f"This can happen when:\n"
                   f"• The server is under heavy load\n"
                   f"• Your sequence is particularly complex\n"
                   f"• The model requires more processing time\n\n"
                   f"💡 Suggestions:\n"
                   f"• Try again in a few minutes\n"
                   f"• Use a shorter sequence\n"
                   f"• Switch to OpenFold2 (usually faster)\n"
                   f"• Enable Demo Mode to test the interface"
    }

def poll_alphafold2_result(req_id: str, api_key: str, status_endpoint: str, model_name: str = "AlphaFold2") -> Dict[str, Any]:
    """
    Poll for AlphaFold2 result using NVIDIA Health API specific endpoint
    """
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "NVCF-POLL-SECONDS": "300"
    }
    
    max_attempts = 180  # 30 minutes at 10-second intervals
    progress_bar = st.progress(0)
    status_placeholder = st.empty()
    
    for attempt in range(max_attempts):
        try:
            progress = attempt / max_attempts
            progress_bar.progress(progress)
            
            minutes_elapsed = (attempt * 10) // 60
            seconds_elapsed = (attempt * 10) % 60
            
            status_placeholder.info(
                f"🔬 {model_name} is processing your protein structure...\n"
                f"Time elapsed: {minutes_elapsed}m {seconds_elapsed}s | "
                f"Request ID: {req_id}"
            )
            
            # Poll using the Health API status endpoint
            poll_response = requests.get(
                f"{status_endpoint}/{req_id}",
                headers=headers,
                timeout=30
            )
            
            if poll_response.status_code == 200:
                # Success - result is ready
                progress_bar.progress(1.0)
                status_placeholder.success(f"🎉 {model_name} structure prediction completed!")
                return {"status": "success", "data": poll_response.json()}
                
            elif poll_response.status_code == 202:
                # Still processing
                if attempt < 5:
                    status_placeholder.info(f"⏳ {model_name} request is queued...")
                elif attempt < 30:
                    status_placeholder.info(f"🧬 {model_name} is running MSA search and structure prediction...")
                else:
                    status_placeholder.info(f"🔬 {model_name} is refining the structure prediction...")
                
                time.sleep(10)
                continue
                
            else:
                # Error status
                error_text = poll_response.text
                return {"status": "error", "message": f"{model_name} failed with status {poll_response.status_code}: {error_text}"}
                
        except requests.exceptions.RequestException as e:
            if attempt == max_attempts - 1:
                return {"status": "error", "message": f"{model_name} polling failed after {max_attempts} attempts: {str(e)}"}
            time.sleep(10)
            continue
    
    timeout_minutes = (max_attempts * 10) // 60
    return {
        "status": "error", 
        "message": f"⏱️ {model_name} prediction timed out after {timeout_minutes} minutes.\n\n"
                   f"AlphaFold2 can take a long time for complex sequences due to:\n"
                   f"• Multiple sequence alignment (MSA) search\n"
                   f"• Template search and processing\n"
                   f"• Neural network structure prediction\n"
                   f"• Structure refinement and relaxation\n\n"
                   f"💡 Suggestions:\n"
                   f"• Try again with a shorter sequence\n"
                   f"• Switch to OpenFold2 (usually faster)\n"
                   f"• Use Demo Mode to test the interface"
    }
def extract_pdb_from_response(response_data: Any) -> Optional[str]:
    """
    Extract PDB content from API response
    """
    if isinstance(response_data, dict):
        # Check for structures_in_ranked_order format (common in protein folding APIs)
        if "structures_in_ranked_order" in response_data:
            structures = response_data["structures_in_ranked_order"]
            if isinstance(structures, list) and len(structures) > 0:
                first_structure = structures[0]
                if isinstance(first_structure, dict) and "structure" in first_structure:
                    pdb_content = first_structure["structure"]
                    if isinstance(pdb_content, str) and ("ATOM" in pdb_content or "HEADER" in pdb_content):
                        return pdb_content
        
        # Try various possible keys where PDB data might be stored
        pdb_keys = ["pdb", "structure", "output", "result", "prediction"]
        
        for key in pdb_keys:
            if key in response_data:
                value = response_data[key]
                if isinstance(value, str) and ("ATOM" in value or "HEADER" in value):
                    return value
                elif isinstance(value, dict):
                    # Recursively search in nested dictionaries
                    nested_pdb = extract_pdb_from_response(value)
                    if nested_pdb:
                        return nested_pdb
                elif isinstance(value, list) and len(value) > 0:
                    # Check if it's a list of structures
                    for item in value:
                        if isinstance(item, dict):
                            nested_pdb = extract_pdb_from_response(item)
                            if nested_pdb:
                                return nested_pdb
        
        # If no PDB found, return the raw response as a string for debugging
        return str(response_data)
    
    elif isinstance(response_data, str):
        if "ATOM" in response_data or "HEADER" in response_data:
            return response_data
    
    return None

def validate_pdb_content(pdb_content: str) -> dict:
    """
    Validate and analyze PDB content
    """
    if not pdb_content or not isinstance(pdb_content, str):
        return {"valid": False, "error": "No PDB content provided"}
    
    lines = pdb_content.split('\n')
    atom_lines = [line for line in lines if line.startswith('ATOM')]
    
    if not atom_lines:
        return {"valid": False, "error": "No ATOM records found in PDB content"}
    
    try:
        # Basic validation - check if atom lines have proper format
        for line in atom_lines[:5]:  # Check first 5 lines
            if len(line) < 54:  # Minimum length for ATOM record
                return {"valid": False, "error": f"Invalid ATOM record format: {line}"}
        
        # Extract basic statistics
        residues = set()
        atoms_count = len(atom_lines)
        
        for line in atom_lines:
            if len(line) > 27:
                residues.add(line[22:27].strip())
        
        return {
            "valid": True,
            "atoms_count": atoms_count,
            "residues_count": len(residues),
            "lines_total": len(lines)
        }
    
    except Exception as e:
        return {"valid": False, "error": f"PDB validation error: {str(e)}"}

# def create_3d_visualization(pdb_content: str) -> str:
#     """
#     Create a 3D molecular visualization using py3Dmol
#     """
#     try:
#         viewer = py3Dmol.view(width=800, height=600)
#         viewer.addModel(pdb_content, 'pdb')
#         viewer.setStyle({'cartoon': {'color': 'spectrum'}})
#         viewer.zoomTo()
#         viewer.spin(True)
#         return viewer._make_html()
#     except Exception as e:
#         st.error(f"Visualization error: {str(e)}")
#         return f"<p>Visualization failed: {str(e)}</p>"

def create_3d_visualization(
    pdb_content: str,
    vmin: float = 0.0,    # lower bound of pLDDT range for color scale
    vmax: float = 100.0,  # upper bound of pLDDT range for color scale
    color_by_plddt: bool = None  # None = auto-detect, True = force pLDDT, False = force solid color
) -> str:
    """
    Create a 3D molecular visualization using py3Dmol, colored by pLDDT
    (stored in the B-factor column, as in AlphaFold outputs).
    Returns HTML string with vertical color bar legend.
    """
    try:
        # Only auto-detect if color_by_plddt is None
        # If explicitly set to True or False, use that value
        if color_by_plddt is None:
            # Check if B-factors look like pLDDT scores
            b_factors = []
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM'):
                    try:
                        b_factor = float(line[60:66].strip())
                        b_factors.append(b_factor)
                    except (ValueError, IndexError):
                        continue
            
            if b_factors:
                min_b = min(b_factors)
                max_b = max(b_factors)
                avg_b = sum(b_factors) / len(b_factors)
                has_variation = (max_b - min_b) > 5
                in_plddt_range = 0 <= min_b <= 100 and 0 <= max_b <= 100
                reasonable_avg = 20 <= avg_b <= 100
                color_by_plddt = has_variation and in_plddt_range and reasonable_avg
            else:
                color_by_plddt = False
        
        # Build custom HTML with 3Dmol viewer and vertical color bar
        import random
        view_id = f"viewer_{random.randint(1000, 9999)}"
        
        # Escape backticks in PDB content for JavaScript
        pdb_escaped = pdb_content.replace('`', '\\`').replace('${', '\\${')
        
        if color_by_plddt:
            style_js = """
                viewer.setStyle({}, {
                    cartoon: {
                        colorfunc: function(atom) {
                            var plddt = atom.b;
                            // AlphaFold-style coloring: purple (low) -> blue -> green -> yellow (high)
                            if (plddt >= 90) {
                                // Yellow to light yellow (90-100)
                                return 'rgb(255, 255, 0)';
                            } else if (plddt >= 70) {
                                // Green to yellow-green (70-90)
                                var t = (plddt - 70) / 20;
                                var r = Math.round(100 + 155 * t);
                                var g = Math.round(200 + 55 * t);
                                var b = Math.round(50 * (1 - t));
                                return 'rgb(' + r + ',' + g + ',' + b + ')';
                            } else if (plddt >= 50) {
                                // Cyan/teal to green (50-70)
                                var t = (plddt - 50) / 20;
                                var r = Math.round(50 + 50 * t);
                                var g = Math.round(180 + 20 * t);
                                var b = Math.round(150 - 100 * t);
                                return 'rgb(' + r + ',' + g + ',' + b + ')';
                            } else {
                                // Purple/blue to cyan (0-50)
                                var t = plddt / 50;
                                var r = Math.round(80 - 30 * t);
                                var g = Math.round(50 + 130 * t);
                                var b = Math.round(150);
                                return 'rgb(' + r + ',' + g + ',' + b + ')';
                            }
                        }
                    }
                });
            """
            legend_html = """
            <div style="
                position: absolute;
                right: 20px;
                top: 50%;
                transform: translateY(-50%);
                display: flex;
                flex-direction: row;
                align-items: center;
                gap: 8px;
                font-family: Arial, sans-serif;
            ">
                <div style="
                    width: 25px;
                    height: 300px;
                    background: linear-gradient(to bottom, 
                        #FFFF00 0%,
                        #90EE90 25%,
                        #32CD9A 50%,
                        #4682B4 75%,
                        #483D8B 100%
                    );
                    border-radius: 4px;
                    border: 1px solid #ccc;
                "></div>
                <div style="
                    display: flex;
                    flex-direction: column;
                    justify-content: space-between;
                    height: 300px;
                    font-size: 12px;
                    color: #333;
                ">
                    <span style="font-weight: bold;">100</span>
                    <span>90</span>
                    <span>70</span>
                    <span>50</span>
                    <span style="font-weight: bold;">0</span>
                </div>
                <div style="
                    writing-mode: vertical-rl;
                    text-orientation: mixed;
                    transform: rotate(180deg);
                    font-size: 11px;
                    color: #555;
                    letter-spacing: 1px;
                ">Prediction Score (pLDDT)</div>
            </div>
            """
        else:
            # Use a single solid color (steel blue) for structures without pLDDT
            style_js = """
                viewer.setStyle({}, {
                    cartoon: { color: '#4682B4' }
                });
            """
            # No legend needed for single color
            legend_html = """
            <div style="
                position: absolute;
                right: 20px;
                bottom: 20px;
                background: rgba(255,255,255,0.9);
                padding: 10px 15px;
                border-radius: 6px;
                font-family: Arial, sans-serif;
                font-size: 12px;
                color: #333;
                border: 1px solid #ccc;
            ">
                <span style="display: inline-block; width: 12px; height: 12px; background: #4682B4; border-radius: 2px; margin-right: 8px; vertical-align: middle;"></span>
                Scaffold Structure
            </div>
            """
        
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ background: #fff; overflow: hidden; }}
        .container {{ position: relative; width: 100%; height: 600px; }}
        #{view_id} {{ width: 100%; height: 100%; }}
        .atom-label {{
            background: rgba(0,0,0,0.8);
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-family: monospace;
            font-size: 11px;
            pointer-events: none;
        }}
        .controls {{
            position: absolute;
            top: 10px;
            left: 10px;
            display: flex;
            flex-direction: column;
            gap: 5px;
            z-index: 100;
        }}
        .ctrl-btn {{
            background: rgba(255,255,255,0.9);
            border: 1px solid #ccc;
            border-radius: 4px;
            padding: 6px 10px;
            cursor: pointer;
            font-size: 12px;
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .ctrl-btn:hover {{
            background: #f0f0f0;
        }}
        .ctrl-btn.active {{
            background: #76B900;
            color: white;
            border-color: #76B900;
        }}
        .info-panel {{
            position: absolute;
            bottom: 10px;
            left: 10px;
            background: rgba(0,0,0,0.85);
            color: white;
            padding: 10px 15px;
            border-radius: 8px;
            font-family: 'Courier New', monospace;
            font-size: 12px;
            max-width: 300px;
            display: none;
            z-index: 100;
        }}
        .info-panel.visible {{
            display: block;
        }}
        .info-title {{
            color: #76B900;
            font-weight: bold;
            margin-bottom: 5px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="controls">
            <button class="ctrl-btn" id="spin-btn" title="Toggle Spin">
                🔄 <span id="spin-text">Spin: ON</span>
            </button>
            <button class="ctrl-btn" id="style-btn" title="Toggle Style">
                🎨 <span id="style-text">Cartoon</span>
            </button>
            <button class="ctrl-btn" id="click-btn" title="Click Mode">
                👆 <span id="click-text">Click: Inspect</span>
            </button>
            <button class="ctrl-btn" id="reset-btn" title="Reset View">
                🔍 <span>Reset View</span>
            </button>
        </div>
        <div id="{view_id}"></div>
        <div class="info-panel" id="info-panel">
            <div class="info-title">Residue Info</div>
            <div id="info-content">Click on a residue to see details</div>
        </div>
        {legend_html}
    </div>
    <script>
        $(document).ready(function() {{
            let viewer = $3Dmol.createViewer("{view_id}", {{
                backgroundColor: 'white'
            }});
            
            let pdbData = `{pdb_escaped}`;
            let model = viewer.addModel(pdbData, "pdb");
            
            // State variables
            var spinning = true;
            var currentStyle = 'cartoon';
            var selectedResidue = null;
            var clickMode = 'inspect';
            
            // Amino acid 3-letter to full name
            var aaNames = {{
                'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine', 'ASP': 'Aspartic Acid',
                'CYS': 'Cysteine', 'GLU': 'Glutamic Acid', 'GLN': 'Glutamine', 'GLY': 'Glycine',
                'HIS': 'Histidine', 'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
                'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline', 'SER': 'Serine',
                'THR': 'Threonine', 'TRP': 'Tryptophan', 'TYR': 'Tyrosine', 'VAL': 'Valine'
            }};
            
            // Number of neighboring residues to show on each side
            var neighborRadius = 5;
            
            // Custom color scheme: green carbons, red oxygens, blue nitrogens, yellow sulfur
            var elementColors = {{
                'C': '#20B2AA',   // Light sea green / teal for carbons
                'N': '#0000FF',   // Blue for nitrogens
                'O': '#FF0000',   // Red for oxygens
                'S': '#FFFF00',   // Yellow for sulfur
                'H': '#FFFFFF',   // White for hydrogens
                'P': '#FFA500'    // Orange for phosphorus
            }};
            
            // Apply initial style
            function applyStyle() {{
                viewer.setStyle({{}}, {{}});
                
                if (currentStyle === 'cartoon') {{
                    {style_js}
                }} else if (currentStyle === 'stick') {{
                    viewer.setStyle({{}}, {{
                        stick: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            radius: 0.15 
                        }},
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            scale: 0.25 
                        }}
                    }});
                }} else if (currentStyle === 'sphere') {{
                    viewer.setStyle({{}}, {{
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }}
                        }}
                    }});
                }}
                
                // Re-highlight selected residue region if any
                if (selectedResidue) {{
                    highlightRegion(selectedResidue.chain, selectedResidue.resi);
                }}
                
                viewer.render();
            }}
            
            // Find and display hydrogen bonds in a region
            function showHydrogenBonds(chain, centerResi, radius) {{
                var startResi = Math.max(1, centerResi - radius);
                var endResi = centerResi + radius;
                
                // Get all atoms in the region - use a broader selection
                var atoms = model.selectedAtoms({{}});
                var regionAtoms = atoms.filter(function(a) {{
                    return a.chain === chain && a.resi >= startResi && a.resi <= endResi;
                }});
                
                // Find potential H-bond donors (N) and acceptors (O)
                var donors = [];
                var acceptors = [];
                
                regionAtoms.forEach(function(atom) {{
                    // Backbone N is donor, backbone O is acceptor
                    if (atom.atom === 'N') {{
                        donors.push(atom);
                    }}
                    if (atom.atom === 'O') {{
                        acceptors.push(atom);
                    }}
                }});
                
                // Check for H-bonds (distance between 2.5 and 3.5 Å)
                donors.forEach(function(donor) {{
                    acceptors.forEach(function(acceptor) {{
                        // Don't bond to self or adjacent residues (i, i+1)
                        var resiDiff = Math.abs(donor.resi - acceptor.resi);
                        if (resiDiff < 2) return;
                        
                        var dx = donor.x - acceptor.x;
                        var dy = donor.y - acceptor.y;
                        var dz = donor.z - acceptor.z;
                        var dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                        
                        // H-bond distance typically 2.5-3.5 Å
                        if (dist >= 2.5 && dist <= 3.5) {{
                            // Create dashed line by drawing multiple small cylinders
                            var numDashes = 8;
                            var dashFraction = 0.6;  // portion of segment that is visible
                            
                            for (var i = 0; i < numDashes; i++) {{
                                var t1 = i / numDashes;
                                var t2 = (i + dashFraction) / numDashes;
                                
                                var x1 = donor.x + t1 * (acceptor.x - donor.x);
                                var y1 = donor.y + t1 * (acceptor.y - donor.y);
                                var z1 = donor.z + t1 * (acceptor.z - donor.z);
                                
                                var x2 = donor.x + t2 * (acceptor.x - donor.x);
                                var y2 = donor.y + t2 * (acceptor.y - donor.y);
                                var z2 = donor.z + t2 * (acceptor.z - donor.z);
                                
                                viewer.addCylinder({{
                                    start: {{x: x1, y: y1, z: z1}},
                                    end: {{x: x2, y: y2, z: z2}},
                                    radius: 0.04,
                                    color: '#00BFFF',
                                    fromCap: 1,
                                    toCap: 1
                                }});
                            }}
                        }}
                    }});
                }});
            }}
            
            // Highlight a region around the clicked residue (±neighborRadius residues)
            function highlightRegion(chain, centerResi) {{
                var startResi = Math.max(1, centerResi - neighborRadius);
                var endResi = centerResi + neighborRadius;
                
                // Add ball-and-stick style for all residues in the region with element colors
                for (var r = startResi; r <= endResi; r++) {{
                    viewer.addStyle({{chain: chain, resi: r}}, {{
                        stick: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            radius: 0.15 
                        }},
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            scale: 0.25 
                        }}
                    }});
                }}
                
                // Highlight the clicked residue with a stronger style and pink outline effect
                viewer.addStyle({{chain: chain, resi: centerResi}}, {{
                    stick: {{ 
                        colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                        radius: 0.25 
                    }},
                    sphere: {{ 
                        colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                        scale: 0.35 
                    }}
                }});
                
                // Show hydrogen bonds in the region
                showHydrogenBonds(chain, centerResi, neighborRadius);
                
                // Add label only to the center (clicked) residue
                viewer.addLabel(chain + centerResi, {{
                    position: {{x: 0, y: 0, z: 0}},
                    backgroundOpacity: 0.8,
                    backgroundColor: 'black',
                    fontColor: 'white',
                    fontSize: 12
                }}, {{chain: chain, resi: centerResi, atom: 'CA'}});
                
                // Zoom camera to the selected region
                viewer.zoomTo({{chain: chain, resi: [startResi, endResi]}}, 800);  // 800ms animation
                
                viewer.render();
            }}
            
            // Clear highlight - reset to base style without any ball-and-stick
            function clearHighlight() {{
                viewer.removeAllShapes();  // Remove H-bond lines
                viewer.removeAllLabels();
                
                // Completely reset all styles to base
                viewer.setStyle({{}}, {{}});
                
                if (currentStyle === 'cartoon') {{
                    {style_js}
                }} else if (currentStyle === 'stick') {{
                    viewer.setStyle({{}}, {{
                        stick: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            radius: 0.15 
                        }},
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            scale: 0.25 
                        }}
                    }});
                }} else if (currentStyle === 'sphere') {{
                    viewer.setStyle({{}}, {{
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }}
                        }}
                    }});
                }}
                
                viewer.render();
            }}
            
            applyStyle();
            viewer.zoomTo();
            viewer.render();
            
            // Spin animation
            function spin() {{
                if (spinning) {{
                    viewer.rotate(0.5, {{x: 0, y: 1, z: 0}});
                    viewer.render();
                }}
                requestAnimationFrame(spin);
            }}
            spin();
            
            // Spin toggle
            $('#spin-btn').click(function() {{
                spinning = !spinning;
                $(this).toggleClass('active', spinning);
                $('#spin-text').text(spinning ? 'Spin: ON' : 'Spin: OFF');
            }});
            
            // Style toggle
            $('#style-btn').click(function() {{
                var styles = ['cartoon', 'stick', 'sphere'];
                var idx = styles.indexOf(currentStyle);
                currentStyle = styles[(idx + 1) % styles.length];
                $('#style-text').text(currentStyle.charAt(0).toUpperCase() + currentStyle.slice(1));
                selectedResidue = null;
                clearHighlight();
                applyStyle();
            }});
            
            // Click mode toggle
            $('#click-btn').click(function() {{
                $(this).toggleClass('active');
                if ($(this).hasClass('active')) {{
                    $('#click-text').text('Click: Active');
                    $('#info-panel').addClass('visible');
                }} else {{
                    $('#click-text').text('Click: Inspect');
                    $('#info-panel').removeClass('visible');
                    // Clear selection FIRST, then clear highlight
                    selectedResidue = null;
                    viewer.removeAllShapes();  // Remove H-bond lines
                    viewer.removeAllLabels();
                    // Re-apply base style without any highlight
                    viewer.setStyle({{}}, {{}});
                    if (currentStyle === 'cartoon') {{
                        {style_js}
                    }} else if (currentStyle === 'stick') {{
                        viewer.setStyle({{}}, {{
                            stick: {{ 
                                colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                                radius: 0.15 
                            }},
                            sphere: {{ 
                                colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                                scale: 0.25 
                            }}
                        }});
                    }} else if (currentStyle === 'sphere') {{
                        viewer.setStyle({{}}, {{
                            sphere: {{ 
                                colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }}
                            }}
                        }});
                    }}
                    // Reset view to show full structure
                    viewer.zoomTo({{}}, 500);
                    viewer.render();
                }}
            }});
            
            // Reset view button
            $('#reset-btn').click(function() {{
                selectedResidue = null;
                viewer.removeAllShapes();
                viewer.removeAllLabels();
                viewer.setStyle({{}}, {{}});
                if (currentStyle === 'cartoon') {{
                    {style_js}
                }} else if (currentStyle === 'stick') {{
                    viewer.setStyle({{}}, {{
                        stick: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            radius: 0.15 
                        }},
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }},
                            scale: 0.25 
                        }}
                    }});
                }} else if (currentStyle === 'sphere') {{
                    viewer.setStyle({{}}, {{
                        sphere: {{ 
                            colorfunc: function(atom) {{ return elementColors[atom.elem] || '#808080'; }}
                        }}
                    }});
                }}
                // Deactivate inspect mode if active
                $('#click-btn').removeClass('active');
                $('#click-text').text('Click: Inspect');
                $('#info-panel').removeClass('visible');
                // Zoom to full structure
                viewer.zoomTo({{}}, 500);
                viewer.render();
            }});
            
            // Click on atom to show details
            viewer.setClickable({{}}, true, function(atom) {{
                if (!$('#click-btn').hasClass('active')) {{
                    // Just toggle spin if not in click mode
                    spinning = !spinning;
                    $('#spin-btn').toggleClass('active', spinning);
                    $('#spin-text').text(spinning ? 'Spin: ON' : 'Spin: OFF');
                    return;
                }}
                
                // Clear previous selection
                clearHighlight();
                
                // Store selection
                selectedResidue = {{
                    chain: atom.chain,
                    resi: atom.resi
                }};
                
                // Highlight the region around the clicked residue
                highlightRegion(atom.chain, atom.resi);
                
                // Get pLDDT from B-factor
                var plddt = atom.b ? atom.b.toFixed(1) : 'N/A';
                var aaFullName = aaNames[atom.resn] || atom.resn;
                
                // Update info panel
                var info = '<b>Residue:</b> ' + atom.resn + ' ' + atom.resi + ' (Chain ' + atom.chain + ')<br>';
                info += '<b>Full Name:</b> ' + aaFullName + '<br>';
                info += '<b>Atom:</b> ' + atom.atom + '<br>';
                info += '<b>pLDDT:</b> ' + plddt + '<br>';
                info += '<b>Coords:</b> (' + atom.x.toFixed(2) + ', ' + atom.y.toFixed(2) + ', ' + atom.z.toFixed(2) + ')';
                
                $('#info-content').html(info);
                $('#info-panel').addClass('visible');
            }});
            
            // Track hovered residue to avoid re-highlighting
            var hoveredResidue = null;
            
            // Hover effect - show pink/magenta highlight like NVIDIA NIMs
            viewer.setHoverable({{}}, true, 
                function(atom, viewer, event, container) {{
                    // Skip if this residue is already hovered
                    var resiKey = atom.chain + '_' + atom.resi;
                    if (hoveredResidue === resiKey) return;
                    
                    hoveredResidue = resiKey;
                    
                    // Add label
                    if (!atom.label) {{
                        atom.label = viewer.addLabel(atom.resn + atom.resi, {{
                            position: atom,
                            backgroundColor: 'rgba(0,0,0,0.8)',
                            fontColor: 'white',
                            fontSize: 11,
                            backgroundOpacity: 0.9
                        }});
                    }}
                    
                    // Add pink/magenta highlight effect on the hovered residue
                    // This creates a glowing outline effect
                    viewer.addStyle({{chain: atom.chain, resi: atom.resi}}, {{
                        cartoon: {{
                            color: '#FF69B4',  // Hot pink for cartoon
                            opacity: 1.0
                        }},
                        stick: {{
                            color: '#FF1493',  // Deep pink for sticks
                            radius: 0.25
                        }},
                        sphere: {{
                            color: '#FF1493',
                            scale: 0.3
                        }}
                    }});
                    
                    viewer.render();
                }},
                function(atom, viewer) {{
                    // Remove label on unhover
                    if (atom.label) {{
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }}
                    
                    var resiKey = atom.chain + '_' + atom.resi;
                    if (hoveredResidue === resiKey) {{
                        hoveredResidue = null;
                        
                        // Re-apply the full style to clear the hover highlight
                        applyStyle();
                    }}
                }}
            );
        }});
    </script>
</body>
</html>
"""
        return html

    except Exception as e:
        return f"<p style='color:red;'>Visualization failed: {str(e)}</p>"

def generate_mock_pdb(sequence: str) -> str:
    """
    Generate a mock PDB structure for demonstration purposes
    """
    pdb_lines = [
        "HEADER    MOCK PROTEIN STRUCTURE                     01-JAN-25   MOCK",
        "TITLE     MOCK PROTEIN STRUCTURE PREDICTION",
        f"SEQRES   1 A {len(sequence):4d}  {' '.join(sequence[:13])}",
    ]
    
    # Add some mock atoms (simplified)
    for i, aa in enumerate(sequence[:50]):  # Limit to first 50 residues for demo
        x, y, z = i * 3.8, 0, 0  # Simple linear arrangement
        pdb_lines.extend([
            f"ATOM  {i*4+1:5d}  N   {aa} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N",
            f"ATOM  {i*4+2:5d}  CA  {aa} A{i+1:4d}    {x+1:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C",
            f"ATOM  {i*4+3:5d}  C   {aa} A{i+1:4d}    {x+2:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C",
            f"ATOM  {i*4+4:5d}  O   {aa} A{i+1:4d}    {x+3:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           O",
        ])
    
    pdb_lines.append("END")
    return "\n".join(pdb_lines)

def main():
    """
    Main Streamlit application
    """
    # NVIDIA Branded Header
    col1, col2 = st.columns([1, 5])
    
    with col2:
        st.markdown("""
        <div style="background: linear-gradient(135deg, #1A1A1A 0%, #2D2D2D 100%); 
                    padding: 20px; border-radius: 8px; border-left: 6px solid #76B900;">
            <h1 style="color: #FFFFFF; margin: 0; border: none;">
                🧬 Protein Structure Prediction
            </h1>
            <p style="color: #76B900; font-weight: 600; margin: 10px 0 0 0;">
                Powered by NVIDIA AI Models
            </p>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    st.markdown("""
    <div style="background-color: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #00D4AA;">
        <p style="margin: 0; color: #1A1A1A;">
            Predict protein 3D structures from amino acid sequences using NVIDIA's state-of-the-art 
            Cloud Functions. Select a model, enter your sequence, and generate accurate structure predictions.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # Sidebar configuration
    try:
        st.sidebar.image("image/nvidia.jpg", use_container_width=True)
    except:
        st.sidebar.markdown("# NVIDIA")
    
    st.sidebar.markdown("---")
    st.sidebar.header("⚙️ Configuration")
    
    # Model selection
    selected_model_name = st.sidebar.selectbox(
        "Select Protein Folding Model:",
        options=list(PROTEIN_MODELS.keys()),
        index=0,
        help="Choose which NVIDIA AI model to use for structure prediction"
    )
    
    selected_model = PROTEIN_MODELS[selected_model_name]
    
    # Show model info with performance hints
    model_info = f"**{selected_model_name}**\n\n{selected_model['description']}"
    
    if "alphafold" in selected_model_name.lower():
        model_info += "\n\n⏳ **Processing Time**: 5-10 minutes\n💡 **Tip**: Very high accuracy but slower"
    elif "openfold" in selected_model_name.lower():
        model_info += "\n\n⚡ **Processing Time**: 2-5 minutes\n💡 **Tip**: Good balance of speed and accuracy"
    elif "boltz" in selected_model_name.lower():
        model_info += "\n\n🚀 **Processing Time**: 3-7 minutes\n💡 **Tip**: Latest improvements"
    
    st.sidebar.info(model_info)
    
    # API configuration - load from environment
    default_api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_CLI_API_KEY") or ""
    api_key = st.sidebar.text_input(
        "NVIDIA API Key:",
        value=default_api_key,
        type="password",
        help="Your NVIDIA Cloud Functions API key"
    )
    
    # Demo mode option
    demo_mode = st.sidebar.checkbox(
        "Demo Mode",
        value=False,
        help="Use mock predictions instead of real API calls for testing"
    )
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("📝 Input")
        
        # Amino Acid Reference Chart
        with st.expander("🧬 Amino Acid Reference Chart", expanded=False):
            st.markdown("**The 20 Standard Amino Acids**")
            
            # Display the amino acid chart image
            try:
                from PIL import Image
                import os
                
                # Look for the amino acid chart image in the current directory
                image_path = None
                for filename in ["amino_acid_chart.webp", "amino_acids.png", "amino_acid_reference.png", 
                               "chart.png", "reference.png", "amino.png"]:
                    if os.path.exists(filename):
                        image_path = filename
                        break
                
                if image_path:
                    image = Image.open(image_path)
                    st.image(image, caption="Amino Acid Reference Chart", use_column_width=True)
                    st.info("💡 **Tip**: Use only the single-letter codes (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y) when entering sequences")
                else:
                    # Fallback if image not found
                    st.warning("⚠️ Amino acid reference image not found. Please add your chart image to the project directory.")
                    st.info("💡 **Expected filenames**: amino_acid_chart.png, amino_acids.png, or similar")
                    st.info("💡 **Single-letter codes**: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
            
            except ImportError:
                st.warning("⚠️ PIL (Pillow) not installed. Please install it to display the amino acid chart image.")
                st.info("💡 **Single-letter codes**: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
            except Exception as e:
                st.warning(f"⚠️ Could not load amino acid chart image: {str(e)}")
                st.info("💡 **Single-letter codes**: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
        
        # Example sequences
        examples = {
            "Select an example...": "",
            "Insulin B-chain (30 AA)": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
            "Lysozyme fragment (140 AA)": "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
            "Sample sequence (27 AA)": "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVST",
            "Cytochrome C fragment (50 AA)": "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"[:50]
        }
        
        selected_example = st.selectbox("Choose an example:", list(examples.keys()))
        
        sequence_input = st.text_area(
            "Enter Protein Amino Acid Sequence:",
            value=examples.get(selected_example, ""),
            height=200,
            placeholder="Enter amino acid sequence using single-letter codes (e.g., ACDEFGHIKLMNPQRSTVWY)",
            help="Enter a protein sequence using standard single-letter amino acid codes"
        )
        
        # Sequence info
        if sequence_input:
            clean_seq = re.sub(r'\s+', '', sequence_input.upper())
            st.info(f"Sequence length: {len(clean_seq)} amino acids")
        
        predict_button = st.button("🔬 Predict Structure", type="primary", disabled=not sequence_input)
    
    with col2:
        st.header("📊 Results")
        
        if predict_button:
            if not api_key and not demo_mode:
                st.error("Please provide an API key or enable Demo Mode")
                return
            
            # Validate sequence
            is_valid, result = validate_protein_sequence(sequence_input)
            
            if not is_valid:
                st.error(f"Invalid sequence: {result}")
                return
            
            clean_sequence = result
            st.success(f"✅ Valid protein sequence with {len(clean_sequence)} amino acids")
            
            # Prediction
            with st.spinner(f"Predicting structure using {selected_model_name}..."):
                try:
                    if demo_mode:
                        # Demo mode - generate mock PDB
                        st.info("Demo mode: Generating mock structure...")
                        time.sleep(2)  # Simulate processing time
                        pdb_content = generate_mock_pdb(clean_sequence)
                        prediction_result = {"status": "success", "data": {"pdb": pdb_content}}
                    else:
                        # Real API call
                        prediction_result = call_nvidia_protein_api(
                            clean_sequence, 
                            selected_model["id"], 
                            api_key,
                            selected_model_name
                        )
                    
                    if prediction_result["status"] == "success":
                        # Extract PDB content
                        pdb_content = extract_pdb_from_response(prediction_result["data"])
                        
                        if pdb_content:
                            # Store in session state
                            st.session_state['pdb_content'] = pdb_content
                            st.session_state['sequence'] = clean_sequence
                            st.session_state['model_used'] = selected_model_name
                            
                            st.success("🎉 Structure prediction completed!")
                        else:
                            st.error("No PDB structure found in the response")
                            st.text("Raw response:")
                            st.text(prediction_result["data"])
                    else:
                        st.error(f"Prediction failed: {prediction_result['message']}")
                        
                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")
    
    # Display results if available
    if 'pdb_content' in st.session_state:
        st.header("🧬 3D Structure Visualization")
        
        # Info about the prediction
        col_info1, col_info2, col_info3 = st.columns(3)
        with col_info1:
            st.metric("Sequence Length", f"{len(st.session_state.get('sequence', ''))} AA")
        with col_info2:
            st.metric("Model Used", st.session_state.get('model_used', 'Unknown'))
        with col_info3:
            st.metric("Status", "✅ Complete")
        
        # Visualization and download
        viz_col1, viz_col2 = st.columns([3, 1])
        
        with viz_col1:
            # Validate PDB content first
            pdb_validation = validate_pdb_content(st.session_state['pdb_content'])
            
            if pdb_validation["valid"]:
                st.success(f"✅ Valid PDB structure with {pdb_validation['atoms_count']} atoms and {pdb_validation['residues_count']} residues")
                
                try:
                    html_content = create_3d_visualization(st.session_state['pdb_content'])
                    components.html(html_content, height=600)
                except Exception as e:
                    st.error(f"❌ 3D Visualization error: {str(e)}")
                    st.info("💡 Try downloading the PDB file and opening it in a molecular viewer like PyMOL or ChimeraX")
                    
                    # Show a text preview instead
                    st.subheader("📋 PDB Text Preview")
                    lines = st.session_state['pdb_content'].split('\n')
                    atom_lines = [line for line in lines if line.startswith('ATOM')][:20]  # Show first 20 atoms
                    preview_text = '\n'.join(atom_lines)
                    if len(atom_lines) == 20:
                        preview_text += f"\n... and {pdb_validation['atoms_count'] - 20} more atoms"
                    st.text_area("First 20 ATOM records:", preview_text, height=300)
            else:
                st.error(f"❌ Invalid PDB content: {pdb_validation['error']}")
                st.info("💡 The API response may not contain valid PDB data. Check the raw content below.")
                
                # Show raw content for debugging
                st.subheader("🔍 Raw API Response")
                st.text_area("Raw Response Content:", st.session_state['pdb_content'][:2000], height=300)
                if len(st.session_state['pdb_content']) > 2000:
                    st.info(f"Showing first 2000 characters of {len(st.session_state['pdb_content'])} total characters")
        
        with viz_col2:
            st.subheader("📥 Download")
            
            filename = f"protein_structure_{len(st.session_state.get('sequence', ''))}aa_{st.session_state.get('model_used', 'unknown').lower().replace(' ', '_')}.pdb"
            
            st.download_button(
                label="📄 Download PDB File",
                data=st.session_state['pdb_content'],
                file_name=filename,
                mime="chemical/x-pdb",
                help="Download the predicted structure as a PDB file"
            )
            
            st.markdown("---")
            st.markdown("**💡 Usage Tips:**")
            st.markdown("- Use the mouse to rotate the structure")
            st.markdown("- Scroll to zoom in/out")
            st.markdown("- The structure shows the predicted 3D conformation")
        
        # Raw PDB content (expandable)
        with st.expander("📄 View Raw PDB Content", expanded=False):
            pdb_content = st.session_state['pdb_content']
            
            # Show PDB statistics
            if pdb_content and isinstance(pdb_content, str) and "ATOM" in pdb_content:
                lines = pdb_content.split('\n')
                atom_lines = [line for line in lines if line.startswith('ATOM')]
                residue_count = len(set(line[22:27].strip() for line in atom_lines if len(line) > 27))
                
                col_stat1, col_stat2, col_stat3 = st.columns(3)
                with col_stat1:
                    st.metric("Total Atoms", len(atom_lines))
                with col_stat2:
                    st.metric("Residues", residue_count)
                with col_stat3:
                    st.metric("PDB Lines", len([l for l in lines if l.strip()]))
                
                st.markdown("**PDB File Content:**")
                st.text_area("", pdb_content, height=400, key="pdb_viewer")
                
                # Add download button for easy access
                st.download_button(
                    label="💾 Download PDB (Alternative)",
                    data=pdb_content,
                    file_name=f"protein_structure_{len(st.session_state.get('sequence', ''))}aa.pdb",
                    mime="chemical/x-pdb",
                    help="Alternative download button for the PDB file"
                )
            else:
                st.warning("⚠️ PDB content appears to be malformed or missing ATOM records")
                st.text_area("Raw Response:", pdb_content, height=400)
    
    # NVIDIA Footer
    st.markdown("---")
    
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col2:
        st.markdown("""
        <div style="text-align: center; color: #666666; padding: 20px;">
            <p style="font-size: 14px; margin-bottom: 10px;">
                <strong>Powered by NVIDIA AI</strong><br>
                Protein structure prediction using NVIDIA Cloud Functions
            </p>
            <p style="font-size: 12px; color: #999999;">
                <strong>Available Models:</strong> OpenFold2, AlphaFold2, AlphaFold2 Multimer, Boltz2<br>
                Predictions are computational estimates. Validate experimentally.<br>
                © 2025 NVIDIA Corporation. All rights reserved.
            </p>
        </div>
        """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
