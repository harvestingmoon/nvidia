#!/usr/bin/env python3
"""
Integration test for loading pipeline results into workflow session
"""
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from workflow.workflow_state import WorkflowSession, WorkflowStage, StageStatus

# Import the loading function from binding_workflow_app
# We'll simulate it here since we can't import streamlit code directly
def load_pipeline_results_test(session: WorkflowSession, folder_path: str) -> bool:
    """Test version of load_pipeline_results without streamlit dependencies"""
    import os
    
    try:
        # Check if folder exists
        if not os.path.exists(folder_path):
            print(f"❌ Folder not found: {folder_path}")
            return False
        
        # Extract project name from folder
        folder_name = os.path.basename(folder_path.rstrip('/'))
        if '_AF2_output' in folder_name:
            project_name = folder_name.replace('_AF2_output', '')
            model_used = 'AlphaFold2'
        elif '_OF3_output' in folder_name:
            project_name = folder_name.replace('_OF3_output', '')
            model_used = 'OpenFold3'
        else:
            project_name = folder_name.replace('_output', '')
            model_used = 'Unknown'
        
        print(f"📁 Loading project: {project_name} (Model: {model_used})")
        
        # Load target structure (Step 1)
        target_pdb = os.path.join(folder_path, f"{project_name}_first_structure.pdb")
        if os.path.exists(target_pdb):
            with open(target_pdb, 'r') as f:
                session.target.pdb_content = f.read()
            session.target.structure_predicted = True
            session.target.structure_file_path = target_pdb
            session.target.prediction_model = model_used
            
            # Extract sequence from PDB (simplified)
            print(f"✅ Loaded target structure ({len(session.target.pdb_content)} bytes)")
            session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
        
        # Load binder scaffold (Step 2)
        scaffold_pdb = os.path.join(folder_path, f"{project_name}_RFD_prediction.pdb")
        if os.path.exists(scaffold_pdb):
            with open(scaffold_pdb, 'r') as f:
                session.binder.scaffold_pdb = f.read()
            session.binder.scaffold_file_path = scaffold_pdb
            print(f"✅ Loaded scaffold structure ({len(session.binder.scaffold_pdb)} bytes)")
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
                    # Parse header
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
            print(f"✅ Loaded {len(sequences)} MPNN sequences (scores: {min(scores):.3f} - {max(scores):.3f})")
            session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.COMPLETED)
        
        # Load complex structure (Step 4)
        complex_pdb = os.path.join(folder_path, f"{project_name}_pdb_1_MULTIMER.pdb")
        if os.path.exists(complex_pdb):
            with open(complex_pdb, 'r') as f:
                session.complex.complex_pdb = f.read()
            print(f"✅ Loaded complex structure ({len(session.complex.complex_pdb)} bytes)")
            session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.COMPLETED)
        
        # Load pLDDT scores (Step 5)
        plddt_file = os.path.join(folder_path, "pLDDT_scores.txt")
        if os.path.exists(plddt_file):
            with open(plddt_file, 'r') as f:
                content = f.read().strip()
            # Parse score (format: "filename.pdb  score")
            parts = content.split()
            for part in parts:
                try:
                    session.complex.plddt_score = float(part.strip())
                    print(f"✅ Loaded pLDDT score: {session.complex.plddt_score:.6f}")
                    break
                except ValueError:
                    continue
        
        # Update project name
        session.project_name = project_name
        
        # Advance to appropriate stage
        if session.complex.complex_pdb:
            session.advance_to_stage(WorkflowStage.RESULTS)
            print(f"✅ Advanced to RESULTS stage")
        
        return True
        
    except Exception as e:
        print(f"❌ Error loading pipeline results: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def test_load_pipeline():
    """Test loading pipeline results into a workflow session"""
    print("\n" + "=" * 70)
    print("🧪 Testing Pipeline Results Loading")
    print("=" * 70)
    
    folders = [
        "workflow/5tpn_AF2_output",
        "workflow/5tpn_OF3_output"
    ]
    
    for folder in folders:
        if not os.path.exists(folder):
            print(f"\n⚠️ Skipping {folder} (not found)")
            continue
        
        print(f"{'=' * 70}")
        print(f"📂 Testing: {folder}")
        print("=" * 70)
        
        # Create new session using the class method
        session = WorkflowSession.create_new(project_name="test_load")
        
        # Check initial state
        print(f"\n📊 Initial State:")
        print(f"   Current Stage: {session.current_stage}")
        print(f"   Target predicted: {session.target.structure_predicted}")
        print(f"   Binder scaffold: {session.binder.scaffold_pdb is not None}")
        print(f"   Sequences count: {len(session.binder.mpnn_sequences)}")
        print(f"   Complex available: {session.complex.complex_pdb is not None}")
        
        # Load results
        print(f"\n🔄 Loading results...")
        success = load_pipeline_results_test(session, folder)
        
        if success:
            print(f"\n✅ Loading succeeded!")
            print(f"\n📊 Final State:")
            print(f"   Project Name: {session.project_name}")
            print(f"   Current Stage: {session.current_stage}")
            print(f"   Target predicted: {session.target.structure_predicted}")
            print(f"   Target model: {session.target.prediction_model}")
            print(f"   Binder scaffold: {session.binder.scaffold_pdb is not None}")
            print(f"   Sequences count: {len(session.binder.mpnn_sequences)}")
            print(f"   Complex available: {session.complex.complex_pdb is not None}")
            if session.complex.plddt_score:
                print(f"   pLDDT score: {session.complex.plddt_score:.6f}")
            
            # Check stage statuses
            print(f"\n📋 Stage Statuses:")
            for stage in WorkflowStage:
                status_str = session.stage_statuses.get(stage.value, StageStatus.NOT_STARTED.value)
                status = StageStatus(status_str)
                icon = "✅" if status == StageStatus.COMPLETED else "⏳" if status == StageStatus.IN_PROGRESS else "⭕"
                print(f"   {icon} {stage.name}: {status.name}")
        else:
            print(f"\n❌ Loading failed!")
    
    print("\n" + "=" * 70)
    print("✅ Test complete!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    test_load_pipeline()
