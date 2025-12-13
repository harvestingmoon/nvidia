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

def create_3d_visualization(pdb_content: str, style: str = "cartoon", show_plddt_legend: bool = True) -> str:
    """Create 3D molecular visualization using py3Dmol with pLDDT coloring"""
    view_id = "mol_view"
    
    # Check if this is AlphaFold output (has B-factor/pLDDT scores)
    # AlphaFold stores pLDDT in the B-factor column
    is_alphafold = "ATOM" in pdb_content and any(
        keyword in pdb_content.upper() 
        for keyword in ["ALPHAFOLD", "PREDICTED", "MODEL"]
    )
    
    style_options = {
        "cartoon": {"cartoon": {"color": "spectrum"}},
        "stick": {"stick": {"radius": 0.2}},
        "sphere": {"sphere": {"radius": 1.0}},
        "line": {"line": {"linewidth": 2}},
        "cartoon+stick": [{"cartoon": {"color": "spectrum"}}, {"stick": {"radius": 0.1}}],
        "plddt": {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}}
    }
    
    selected_style = style_options.get(style, style_options["cartoon"])
    
    # Legend HTML for pLDDT scores
    legend_html = ""
    if show_plddt_legend:
        legend_html = """
        <div style="position: absolute; bottom: 10px; right: 10px; background: rgba(255,255,255,0.95); 
                    padding: 12px 15px; border-radius: 8px; font-family: Arial, sans-serif; font-size: 12px;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.15); z-index: 1000; border: 1px solid #e0e0e0;">
            <div style="font-weight: bold; margin-bottom: 8px; color: #333; font-size: 13px;">
                pLDDT Confidence
            </div>
            <div style="display: flex; flex-direction: column; gap: 4px;">
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #0053D6; border-radius: 2px;"></div>
                    <span style="color: #444;">&gt;90: Very high</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #65CBF3; border-radius: 2px;"></div>
                    <span style="color: #444;">70-90: Confident</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #FFDB13; border-radius: 2px;"></div>
                    <span style="color: #444;">50-70: Low</span>
                </div>
                <div style="display: flex; align-items: center; gap: 8px;">
                    <div style="width: 20px; height: 14px; background: #FF7D45; border-radius: 2px;"></div>
                    <span style="color: #444;">&lt;50: Very low</span>
                </div>
            </div>
        </div>
        """
    
    html_content = f"""
    <div id="{view_id}" style="height: 600px; width: 100%; position: relative;">
        {legend_html}
    </div>
    <script src="https://cdn.jsdelivr.net/npm/3dmol@latest/build/3Dmol-min.js"></script>
    <script>
    $(document).ready(function() {{
        let viewer = $3Dmol.createViewer($("#{view_id}"), {{
            defaultcolors: $3Dmol.rasmolElementColors
        }});
        
        let pdbData = `{pdb_content}`;
        
        viewer.addModel(pdbData, "pdb");
        
        // Apply pLDDT coloring using B-factor values (AlphaFold stores pLDDT there)
        viewer.setStyle({{}}, {{
            cartoon: {{
                colorfunc: function(atom) {{
                    // pLDDT is stored in B-factor column
                    var plddt = atom.b;
                    if (plddt > 90) return '#0053D6';      // Dark blue - very high
                    else if (plddt > 70) return '#65CBF3'; // Light blue - confident  
                    else if (plddt > 50) return '#FFDB13'; // Yellow - low
                    else return '#FF7D45';                  // Orange - very low
                }}
            }}
        }});
        
        viewer.zoomTo();
        viewer.render();
        viewer.zoom(1.2);
        
        // Add rotation animation
        viewer.rotate(10, {{x:1, y:1, z:0}});
    }});
    </script>
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
