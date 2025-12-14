#!/usr/bin/env python3
"""
Standalone PDB File Viewer Utility
Simple tool to view and analyze PDB files from protein structure predictions
"""

import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
from typing import Optional, Dict, Any
import os

def validate_pdb_content(pdb_content: str) -> dict:
    """Validate and analyze PDB content"""
    if not pdb_content or not isinstance(pdb_content, str):
        return {"valid": False, "error": "No PDB content provided"}
    
    lines = pdb_content.split('\n')
    atom_lines = [line for line in lines if line.startswith('ATOM')]
    
    if not atom_lines:
        return {"valid": False, "error": "No ATOM records found in PDB content"}
    
    try:
        # Extract detailed statistics
        residues = {}
        atoms_by_residue = {}
        atom_types = set()
        chains = set()
        
        for line in atom_lines:
            if len(line) > 54:
                chain = line[21] if len(line) > 21 else 'A'
                res_num = line[22:26].strip()
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                
                chains.add(chain)
                atom_types.add(atom_name)
                
                res_key = f"{chain}:{res_num}:{res_name}"
                residues[res_key] = res_name
                
                if res_name not in atoms_by_residue:
                    atoms_by_residue[res_name] = 0
                atoms_by_residue[res_name] += 1
        
        return {
            "valid": True,
            "atoms_count": len(atom_lines),
            "residues_count": len(residues),
            "unique_residue_types": len(set(residues.values())),
            "chains": list(chains),
            "atom_types": list(atom_types),
            "residue_composition": atoms_by_residue,
            "lines_total": len(lines)
        }
    
    except Exception as e:
        return {"valid": False, "error": f"PDB validation error: {str(e)}"}


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
    # If all values are the same or outside 0-100, probably not pLDDT
    has_variation = (max_b - min_b) > 5
    in_plddt_range = 0 <= min_b <= 100 and 0 <= max_b <= 100
    reasonable_avg = 20 <= avg_b <= 100
    
    return has_variation and in_plddt_range and reasonable_avg

def create_3d_visualization(pdb_content: str, style: str = "cartoon", show_plddt_legend: bool = True, color_by_plddt: bool = None) -> str:
    """Create 3D molecular visualization using py3Dmol with pLDDT coloring
    
    Args:
        pdb_content: PDB file content as string
        style: Visualization style ('cartoon', 'stick', etc.)
        show_plddt_legend: Whether to show the pLDDT legend
        color_by_plddt: Force pLDDT coloring on/off. If None, auto-detect.
    """
    import random
    view_id = f"mol_view_{random.randint(1000, 9999)}"
    
    # Auto-detect if we should use pLDDT coloring
    if color_by_plddt is None:
        color_by_plddt = check_has_plddt_scores(pdb_content)
    
    # Show legend only if we're using pLDDT coloring
    show_legend = show_plddt_legend and color_by_plddt
    
    # Legend HTML for pLDDT scores - using absolute positioning within container
    legend_html = ""
    if show_legend:
        legend_html = f"""
        <div id="legend_{view_id}" style="
            position: absolute; 
            top: 10px; 
            right: 10px; 
            background: rgba(30, 30, 30, 0.95); 
            padding: 12px 16px; 
            border-radius: 10px; 
            font-family: 'Segoe UI', Arial, sans-serif; 
            font-size: 11px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.4); 
            z-index: 1000; 
            border: 1px solid rgba(118, 185, 0, 0.5);
            pointer-events: none;
        ">
            <div style="font-weight: 600; margin-bottom: 8px; color: #76B900; font-size: 12px; letter-spacing: 0.5px;">
                pLDDT Confidence
            </div>
            <div style="display: flex; flex-direction: column; gap: 5px;">
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #0053D6; border-radius: 2px;"></div>
                    <span style="color: #E8E8E8;">&gt;90: Very high</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #65CBF3; border-radius: 2px;"></div>
                    <span style="color: #E8E8E8;">70-90: Confident</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #FFDB13; border-radius: 2px;"></div>
                    <span style="color: #E8E8E8;">50-70: Low</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #FF7D45; border-radius: 2px;"></div>
                    <span style="color: #E8E8E8;">&lt;50: Very low</span>
                </div>
            </div>
        </div>
        """
    
    # Determine coloring style based on whether pLDDT is available
    if color_by_plddt:
        style_js = """
                viewer.setStyle({}, {
                    cartoon: {
                        colorfunc: function(atom) {
                            var plddt = atom.b;
                            if (plddt > 90) return '#0053D6';
                            else if (plddt > 70) return '#65CBF3';
                            else if (plddt > 50) return '#FFDB13';
                            else return '#FF7D45';
                        }
                    }
                });
        """
    else:
        # Use spectrum coloring by residue number when pLDDT not available
        style_js = """
                viewer.setStyle({}, {
                    cartoon: {
                        color: 'spectrum'
                    }
                });
        """
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/3dmol@latest/build/3Dmol-min.js"></script>
        <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            body {{ margin: 0; padding: 0; overflow: hidden; background: #1a1a1a; }}
            .viewer-container {{ position: relative; width: 100%; height: 500px; }}
            #{view_id} {{ width: 100%; height: 100%; }}
        </style>
    </head>
    <body>
        <div class="viewer-container">
            {legend_html}
            <div id="{view_id}"></div>
        </div>
        <script>
            $(document).ready(function() {{
                let viewer = $3Dmol.createViewer("{view_id}", {{
                    backgroundColor: '#1a1a1a'
                }});
                
                let pdbData = `{pdb_content}`;
                
                viewer.addModel(pdbData, "pdb");
                
                {style_js}
                
                viewer.zoomTo();
                viewer.render();
                viewer.zoom(1.2);
                viewer.rotate(10, {{x:1, y:1, z:0}});
            }});
        </script>
    </body>
    </html>
    """
    
    return html_content

def analyze_sequence_from_pdb(pdb_content: str) -> Optional[str]:
    """Extract amino acid sequence from PDB ATOM records"""
    if not pdb_content:
        return None
    
    lines = pdb_content.split('\n')
    atom_lines = [line for line in lines if line.startswith('ATOM')]
    
    # Standard 3-letter to 1-letter amino acid code mapping
    aa_mapping = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    try:
        residues = {}
        for line in atom_lines:
            if len(line) > 26:
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                if res_name in aa_mapping:
                    residues[res_num] = aa_mapping[res_name]
        
        # Sort by residue number and concatenate
        sorted_residues = sorted(residues.items())
        sequence = ''.join([res[1] for res in sorted_residues])
        
        return sequence
    except:
        return None

def main():
    st.set_page_config(
        page_title="PDB Viewer",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 PDB File Viewer")
    st.markdown("Upload and visualize protein structure files in PDB format")
    
    # Sidebar options
    st.sidebar.header("⚙️ Viewer Options")
    
    visualization_style = st.sidebar.selectbox(
        "Visualization Style:",
        ["cartoon", "stick", "sphere", "line", "cartoon+stick"],
        index=0
    )
    
    # File upload option
    uploaded_file = st.file_uploader("Choose a PDB file", type=['pdb', 'txt'])
    
    # Text input option
    st.subheader("📝 Or Paste PDB Content")
    pdb_text = st.text_area(
        "Paste PDB content here:",
        height=150,
        help="Paste the complete PDB file content including HEADER and ATOM records"
    )
    
    pdb_content = None
    
    # Process uploaded file
    if uploaded_file is not None:
        try:
            pdb_content = str(uploaded_file.read(), "utf-8")
            st.success(f"✅ Loaded PDB file: {uploaded_file.name}")
        except Exception as e:
            st.error(f"❌ Error reading file: {str(e)}")
    
    # Process text input
    elif pdb_text.strip():
        pdb_content = pdb_text.strip()
        st.success("✅ PDB content loaded from text input")
    
    if pdb_content:
        # Validate PDB content
        validation = validate_pdb_content(pdb_content)
        
        if validation["valid"]:
            # Show structure information
            st.subheader("📊 Structure Information")
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Atoms", validation["atoms_count"])
            with col2:
                st.metric("Residues", validation["residues_count"])
            with col3:
                st.metric("Chains", len(validation["chains"]))
            with col4:
                st.metric("Residue Types", validation["unique_residue_types"])
            
            # Extract sequence
            sequence = analyze_sequence_from_pdb(pdb_content)
            if sequence:
                st.subheader("🧬 Amino Acid Sequence")
                st.code(sequence, language="text")
                st.info(f"Sequence Length: {len(sequence)} amino acids")
            
            # 3D Visualization
            st.subheader("🎮 3D Structure Visualization")
            
            try:
                html_content = create_3d_visualization(pdb_content, visualization_style)
                components.html(html_content, height=650)
                
                st.markdown("""
                **💡 Interaction Tips:**
                - **Rotate**: Click and drag
                - **Zoom**: Mouse wheel or pinch
                - **Pan**: Right-click and drag
                """)
                
            except Exception as e:
                st.error(f"❌ Visualization error: {str(e)}")
                st.info("💡 Try a different visualization style or check the PDB format")
            
            # Detailed Analysis
            with st.expander("🔍 Detailed Analysis", expanded=False):
                st.json(validation)
                
                st.subheader("Chain Information")
                for chain in validation["chains"]:
                    st.write(f"**Chain {chain}**")
                
                st.subheader("Residue Composition")
                if validation["residue_composition"]:
                    import pandas as pd
                    df = pd.DataFrame(list(validation["residue_composition"].items()), 
                                    columns=["Residue", "Count"])
                    st.bar_chart(df.set_index("Residue"))
            
            # Raw PDB Content
            with st.expander("📄 Raw PDB Content"):
                st.text_area("PDB File Content:", pdb_content, height=400, key="raw_pdb")
                
                # Download button
                st.download_button(
                    label="💾 Download PDB File",
                    data=pdb_content,
                    file_name="protein_structure.pdb",
                    mime="chemical/x-pdb"
                )
        
        else:
            st.error(f"❌ Invalid PDB content: {validation['error']}")
            
            # Show problematic content for debugging
            st.subheader("🔍 Content Preview")
            st.text_area("First 1000 characters:", pdb_content[:1000], height=200)
    
    else:
        st.info("👆 Please upload a PDB file or paste PDB content to get started")
        
        # Show example
        st.subheader("📋 Example Usage")
        st.markdown("""
        **PDB files typically contain:**
        - `HEADER` records with metadata
        - `ATOM` records with atomic coordinates
        - `TER` and `END` records for structure termination
        
        **Supported formats:**
        - Standard PDB format (.pdb files)
        - Text files containing PDB data
        - Pasted PDB content
        """)
        
        # Example PDB content
        example_pdb = """HEADER    EXAMPLE PROTEIN                         01-JAN-70   1ABC
ATOM      1  N   ALA A   1      20.154   1.850  22.430  1.00 10.00           N
ATOM      2  CA  ALA A   1      19.030   2.650  22.946  1.00 10.00           C
ATOM      3  C   ALA A   1      17.774   1.859  23.296  1.00 10.00           C
ATOM      4  O   ALA A   1      17.795   0.780  23.877  1.00 10.00           O
END"""
        
        if st.button("📥 Load Example PDB"):
            st.rerun()

if __name__ == "__main__":
    main()
