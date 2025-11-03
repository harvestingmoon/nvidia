#!/usr/bin/env python3
"""
Test script to verify folder loading functionality
"""
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from workflow.workflow_state import WorkflowSession, WorkflowStage

def test_folder_structure(folder_path):
    """Test if folder has expected structure"""
    print(f"\n🔍 Testing folder: {folder_path}")
    print("=" * 60)
    
    if not os.path.exists(folder_path):
        print(f"❌ Folder does not exist: {folder_path}")
        return False
    
    # Extract project name
    folder_name = os.path.basename(folder_path.rstrip('/'))
    if '_AF2_output' in folder_name:
        project_name = folder_name.replace('_AF2_output', '')
        model = 'AlphaFold2'
    elif '_OF3_output' in folder_name:
        project_name = folder_name.replace('_OF3_output', '')
        model = 'OpenFold3'
    else:
        project_name = folder_name.replace('_output', '')
        model = 'Unknown'
    
    print(f"📁 Project name: {project_name}")
    print(f"🤖 Model: {model}")
    
    # Check for expected files
    expected_files = {
        'target_structure': f"{project_name}_first_structure.pdb",
        'scaffold': f"{project_name}_RFD_prediction.pdb",
        'sequences': f"{project_name}_Protein_MPNN_prediction.fa",
        'complex': f"{project_name}_pdb_1_MULTIMER.pdb",
        'scores': "pLDDT_scores.txt"
    }
    
    files_found = {}
    for key, filename in expected_files.items():
        filepath = os.path.join(folder_path, filename)
        exists = os.path.exists(filepath)
        files_found[key] = exists
        
        icon = "✅" if exists else "❌"
        size = os.path.getsize(filepath) if exists else 0
        print(f"{icon} {key:20} {filename:50} ({size:,} bytes)")
    
    all_found = all(files_found.values())
    print(f"\n{'✅ All files found!' if all_found else '⚠️ Some files missing'}")
    
    return all_found

def test_file_parsing(folder_path):
    """Test parsing content from files"""
    print(f"\n📖 Testing file parsing")
    print("=" * 60)
    
    folder_name = os.path.basename(folder_path.rstrip('/'))
    project_name = folder_name.replace('_AF2_output', '').replace('_OF3_output', '')
    
    # Test target PDB parsing
    target_pdb = os.path.join(folder_path, f"{project_name}_first_structure.pdb")
    if os.path.exists(target_pdb):
        with open(target_pdb, 'r') as f:
            lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        print(f"✅ Target PDB: {len(lines)} total lines, {len(atom_lines)} ATOM records")
    
    # Test scaffold PDB parsing
    scaffold_pdb = os.path.join(folder_path, f"{project_name}_RFD_prediction.pdb")
    if os.path.exists(scaffold_pdb):
        with open(scaffold_pdb, 'r') as f:
            lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        print(f"✅ Scaffold PDB: {len(lines)} total lines, {len(atom_lines)} ATOM records")
    
    # Test MPNN FASTA parsing
    mpnn_fasta = os.path.join(folder_path, f"{project_name}_Protein_MPNN_prediction.fa")
    if os.path.exists(mpnn_fasta):
        with open(mpnn_fasta, 'r') as f:
            content = f.read()
        sequences = []
        scores = []
        for line in content.split('\n'):
            if line.startswith('>'):
                if 'score=' in line:
                    score_part = line.split('score=')[1].split(',')[0]
                    scores.append(float(score_part))
            elif line.strip() and not line.startswith('>'):
                sequences.append(line.strip())
        
        print(f"✅ MPNN FASTA: {len(sequences)} sequences")
        if sequences:
            print(f"   - First sequence length: {len(sequences[0])} AA")
            if scores:
                print(f"   - Score range: {min(scores):.3f} - {max(scores):.3f}")
    
    # Test complex PDB parsing
    complex_pdb = os.path.join(folder_path, f"{project_name}_pdb_1_MULTIMER.pdb")
    if os.path.exists(complex_pdb):
        with open(complex_pdb, 'r') as f:
            lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        chains = set()
        for line in atom_lines:
            if len(line) > 21:
                chains.add(line[21])
        print(f"✅ Complex PDB: {len(lines)} total lines, {len(atom_lines)} ATOM records, {len(chains)} chains")
    
    # Test pLDDT scores
    plddt_file = os.path.join(folder_path, "pLDDT_scores.txt")
    if os.path.exists(plddt_file):
        with open(plddt_file, 'r') as f:
            content = f.read().strip()
        try:
            score = float(content)
            print(f"✅ pLDDT score: {score:.6f}")
        except ValueError:
            print(f"⚠️ Could not parse pLDDT score: {content}")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("🧪 Pipeline Folder Loading Test")
    print("=" * 60)
    
    # Test both output folders
    folders = [
        "workflow/5tpn_AF2_output",
        "workflow/5tpn_OF3_output"
    ]
    
    for folder in folders:
        if os.path.exists(folder):
            test_folder_structure(folder)
            test_file_parsing(folder)
        else:
            print(f"\n⚠️ Folder not found: {folder}")
    
    print("\n" + "=" * 60)
    print("✅ Test complete!")
    print("=" * 60 + "\n")
