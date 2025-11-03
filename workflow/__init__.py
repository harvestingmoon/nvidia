"""
Workflow Management Package
Contains workflow state management, generative pipeline, and binding analysis modules
"""

from .workflow_state import (
    WorkflowSession,
    WorkflowStage,
    StageStatus,
    WorkflowValidator,
    TargetProteinData,
    BinderProteinData,
    ComplexAnalysisData
)

from .generative_pipeline import GenerativePipeline

from .binding_analysis import (
    parse_pdb_content,
    find_interface_residues,
    assess_binding_quality,
    combine_pdbs,
    generate_contact_map_data
)

__all__ = [
    'WorkflowSession',
    'WorkflowStage',
    'StageStatus',
    'WorkflowValidator',
    'TargetProteinData',
    'BinderProteinData',
    'ComplexAnalysisData',
    'GenerativePipeline',
    'parse_pdb_content',
    'find_interface_residues',
    'assess_binding_quality',
    'combine_pdbs',
    'generate_contact_map_data'
]

__version__ = "1.0.0"
