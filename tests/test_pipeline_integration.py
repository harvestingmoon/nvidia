#!/usr/bin/env python3
"""
Test script for Generative Pipeline integration with WorkflowState
Demonstrates the complete workflow from target input to complex prediction
"""

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from workflow import (
    WorkflowSession,
    WorkflowStage,
    StageStatus,
    GenerativePipeline
)


def test_pipeline_integration():
    """Test the generative pipeline integration"""
    
    print("\n" + "="*70)
    print("GENERATIVE PIPELINE - WORKFLOW STATE INTEGRATION TEST")
    print("="*70 + "\n")
    
    # Get API key from environment (loaded from .env)
    api_key = os.getenv("NGC_CLI_API_KEY") or os.getenv("NVIDIA_API_KEY")
    
    if not api_key:
        print("❌ ERROR: No API key found!")
        print("   Please set NVIDIA_API_KEY in your .env file")
        print("   Copy .env.example to .env and add your key\n")
        return False
    
    print(f"✅ Using API key from environment: {api_key[:20]}...\n")
    
    # Create a new workflow session
    session = WorkflowSession.create_new(project_name="test_insulin")
    print(f"✅ Created workflow session: {session.session_id}")
    print(f"   Project: {session.project_name}")
    print(f"   Current stage: {session.current_stage.value}\n")
    
    # Set target sequence (Insulin B-chain - small for fast testing)
    target_sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    session.target.sequence = target_sequence
    session.target.input_type = "sequence"
    
    print(f"✅ Target sequence set: {len(target_sequence)} amino acids")
    print(f"   Sequence: {target_sequence[:30]}...\n")
    
    # Initialize pipeline
    pipeline = GenerativePipeline(
        session=session,
        api_key=api_key,
        output_dir=Path("test_output")
    )
    
    print(f"✅ Pipeline initialized")
    print(f"   Output directory: {pipeline.output_dir}\n")
    
    # Test stage validation
    print("Testing workflow stage validation:")
    for stage in WorkflowStage:
        can_advance, msg = session.can_advance_to(stage)
        status = "✅" if can_advance else "❌"
        print(f"   {status} {stage.value}: {msg}")
    
    print("\n" + "-"*70)
    print("WORKFLOW STATE SUMMARY")
    print("-"*70)
    summary = pipeline.get_stage_summary()
    for key, value in summary.items():
        print(f"   {key}: {value}")
    
    print("\n" + "-"*70)
    print("SESSION DATA STRUCTURE")
    print("-"*70)
    print(f"   Target:")
    print(f"      - sequence: {bool(session.target.sequence)}")
    print(f"      - pdb_content: {bool(session.target.pdb_content)}")
    print(f"      - model_used: {session.target.model_used}")
    print(f"      - structure_predicted: {session.target.structure_predicted}")
    print(f"   Binder:")
    print(f"      - scaffold_pdb: {bool(session.binder.scaffold_pdb)}")
    print(f"      - mpnn_sequences: {len(session.binder.mpnn_sequences)}")
    print(f"      - design_method: {session.binder.design_method}")
    print(f"   Complex:")
    print(f"      - complex_pdb: {bool(session.complex.complex_pdb)}")
    print(f"      - plddt_score: {session.complex.plddt_score}")
    print(f"      - docking_method: {session.complex.docking_method}")
    
    print("\n" + "="*70)
    print("INTEGRATION TEST COMPLETE")
    print("="*70)
    print("\nTo run the full pipeline with API calls:")
    print("   python test_pipeline_integration.py --run-full")
    print("\nOr run individual steps:")
    print("   python test_pipeline_integration.py --step target")
    print("   python test_pipeline_integration.py --step scaffold")
    print("   python test_pipeline_integration.py --step sequence")
    print("   python test_pipeline_integration.py --step complex")
    
    return True


def run_full_pipeline():
    """Run the complete pipeline (makes actual API calls)"""
    
    api_key = os.getenv("NGC_CLI_API_KEY") or os.getenv("NVIDIA_API_KEY")
    
    if not api_key:
        print("❌ ERROR: No API key found! Set NVIDIA_API_KEY in .env")
        return False
    
    # Create session
    session = WorkflowSession.create_new(project_name="insulin_binder")
    session.target.sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    
    # Initialize pipeline
    pipeline = GenerativePipeline(session, api_key, Path("insulin_binder_output"))
    
    # Run full pipeline
    success, msg = pipeline.run_full_pipeline(
        model="AF2",
        contigs="A1-25/0 70-100",
        hotspot_res=["A14", "A15", "A17", "A18"],
        num_sequences=5,
        num_candidates=1
    )
    
    if success:
        print(f"\n✅ {msg}")
        print(f"\nResults saved to: {pipeline.output_dir}")
        
        # Save session
        session_file = pipeline.output_dir / "session.json"
        with open(session_file, 'w') as f:
            f.write(session.to_json())
        print(f"Session saved to: {session_file}")
    else:
        print(f"\n❌ {msg}")
    
    return success


def run_single_step(step: str):
    """Run a single pipeline step"""
    
    api_key = os.getenv("NGC_CLI_API_KEY") or os.getenv("NVIDIA_API_KEY")
    
    if not api_key:
        print("❌ ERROR: No API key found! Set NVIDIA_API_KEY in .env")
        return False
    
    session = WorkflowSession.create_new(project_name="test_step")
    session.target.sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    
    pipeline = GenerativePipeline(session, api_key, Path(f"test_{step}_output"))
    
    if step == "target":
        success, msg = pipeline.run_target_prediction("AF2")
    elif step == "scaffold":
        # Need to load previous results or mock them
        print("⚠️  This step requires target prediction results")
        print("    Run full pipeline or provide existing PDB")
        return False
    elif step == "sequence":
        print("⚠️  This step requires scaffold design results")
        return False
    elif step == "complex":
        print("⚠️  This step requires sequence design results")
        return False
    else:
        print(f"❌ Unknown step: {step}")
        return False
    
    print(f"\n{'✅' if success else '❌'} {msg}")
    return success


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test Generative Pipeline Integration")
    parser.add_argument("--run-full", action="store_true", help="Run full pipeline with API calls")
    parser.add_argument("--step", choices=["target", "scaffold", "sequence", "complex"],
                       help="Run a single pipeline step")
    
    args = parser.parse_args()
    
    if args.run_full:
        success = run_full_pipeline()
    elif args.step:
        success = run_single_step(args.step)
    else:
        success = test_pipeline_integration()
    
    sys.exit(0 if success else 1)
