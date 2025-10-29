"""
Workflow Management Package
Contains workflow state management and binding analysis modules
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
    'parse_pdb_content',
    'find_interface_residues',
    'assess_binding_quality',
    'combine_pdbs',
    'generate_contact_map_data'
]

__version__ = "1.0.0"
