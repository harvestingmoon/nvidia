"""
Workflow State Management for Binding Protein Design Pipeline
Handles state persistence, validation, and transitions between workflow stages
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, Any, List
from enum import Enum
import json
import time
from datetime import datetime


class WorkflowStage(Enum):
    """Enumeration of workflow stages"""
    TARGET_INPUT = "target_input"
    TARGET_PREDICTION = "target_prediction"
    BINDER_DESIGN = "binder_design"
    BINDER_PREDICTION = "binder_prediction"
    COMPLEX_ANALYSIS = "complex_analysis"
    RESULTS = "results"


class StageStatus(Enum):
    """Status of each stage"""
    NOT_STARTED = "not_started"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class TargetProteinData:
    """Data structure for target protein"""
    sequence: Optional[str] = None
    pdb_content: Optional[str] = None
    pdb_id: Optional[str] = None
    input_type: str = "sequence"  # sequence, pdb_file, pdb_id
    structure_predicted: bool = False
    model_used: Optional[str] = None
    plddt_scores: Optional[List[float]] = None
    confidence_avg: Optional[float] = None
    binding_site_residues: List[int] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BinderProteinData:
    """Data structure for binder protein"""
    sequence: Optional[str] = None
    pdb_content: Optional[str] = None
    design_method: str = "manual"  # manual, rfdiffusion, template
    structure_predicted: bool = False
    model_used: Optional[str] = None
    plddt_scores: Optional[List[float]] = None
    confidence_avg: Optional[float] = None
    mpnn_sequences: List[str] = field(default_factory=list)
    selected_sequence_idx: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ComplexAnalysisData:
    """Data structure for complex analysis results"""
    docking_method: str = "overlay"  # overlay, diffdock
    complex_pdb: Optional[str] = None
    interface_residues_target: List[int] = field(default_factory=list)
    interface_residues_binder: List[int] = field(default_factory=list)
    num_contacts: int = 0
    avg_distance: float = 0.0
    min_distance: float = 0.0
    binding_affinity: Optional[float] = None
    quality_score: int = 0
    quality_grade: str = "N/A"
    feedback: List[str] = field(default_factory=list)
    docking_poses: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class WorkflowSession:
    """Complete workflow session state"""
    session_id: str
    project_name: str
    created_at: str
    last_updated: str
    current_stage: WorkflowStage = WorkflowStage.TARGET_INPUT
    stage_statuses: Dict[str, str] = field(default_factory=dict)
    target: TargetProteinData = field(default_factory=TargetProteinData)
    binder: BinderProteinData = field(default_factory=BinderProteinData)
    complex: ComplexAnalysisData = field(default_factory=ComplexAnalysisData)
    notes: str = ""
    
    def __post_init__(self):
        """Initialize stage statuses"""
        if not self.stage_statuses:
            self.stage_statuses = {
                stage.value: StageStatus.NOT_STARTED.value 
                for stage in WorkflowStage
            }
    
    def update_stage_status(self, stage: WorkflowStage, status: StageStatus):
        """Update the status of a specific stage"""
        self.stage_statuses[stage.value] = status.value
        self.last_updated = datetime.now().isoformat()
    
    def advance_to_stage(self, stage: WorkflowStage):
        """Move to the next stage"""
        self.current_stage = stage
        self.update_stage_status(stage, StageStatus.IN_PROGRESS)
    
    def can_advance_to(self, stage: WorkflowStage) -> tuple[bool, str]:
        """Check if workflow can advance to a specific stage"""
        if stage == WorkflowStage.TARGET_INPUT:
            return True, "Starting workflow"
        
        elif stage == WorkflowStage.TARGET_PREDICTION:
            if self.target.sequence or self.target.pdb_content or self.target.pdb_id:
                return True, "Target input provided"
            return False, "No target protein input provided"
        
        elif stage == WorkflowStage.BINDER_DESIGN:
            if self.target.structure_predicted or self.target.pdb_content:
                return True, "Target structure available"
            return False, "Target structure not available"
        
        elif stage == WorkflowStage.BINDER_PREDICTION:
            if self.binder.sequence:
                return True, "Binder sequence provided"
            return False, "No binder sequence provided"
        
        elif stage == WorkflowStage.COMPLEX_ANALYSIS:
            if (self.target.pdb_content and self.binder.pdb_content):
                return True, "Both structures available"
            return False, "Both target and binder structures required"
        
        elif stage == WorkflowStage.RESULTS:
            if self.complex.complex_pdb or self.complex.num_contacts > 0:
                return True, "Analysis complete"
            return False, "Complex analysis not completed"
        
        return False, "Unknown stage"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert session to dictionary for serialization"""
        return {
            "session_id": self.session_id,
            "project_name": self.project_name,
            "created_at": self.created_at,
            "last_updated": self.last_updated,
            "current_stage": self.current_stage.value,
            "stage_statuses": self.stage_statuses,
            "target": self.target.to_dict(),
            "binder": self.binder.to_dict(),
            "complex": self.complex.to_dict(),
            "notes": self.notes
        }
    
    def to_json(self) -> str:
        """Serialize to JSON string"""
        return json.dumps(self.to_dict(), indent=2)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'WorkflowSession':
        """Create session from dictionary"""
        session = cls(
            session_id=data["session_id"],
            project_name=data["project_name"],
            created_at=data["created_at"],
            last_updated=data["last_updated"],
            current_stage=WorkflowStage(data["current_stage"]),
            stage_statuses=data["stage_statuses"],
            notes=data.get("notes", "")
        )
        
        # Reconstruct nested objects
        session.target = TargetProteinData(**data["target"])
        session.binder = BinderProteinData(**data["binder"])
        session.complex = ComplexAnalysisData(**data["complex"])
        
        return session
    
    @classmethod
    def from_json(cls, json_str: str) -> 'WorkflowSession':
        """Deserialize from JSON string"""
        data = json.loads(json_str)
        return cls.from_dict(data)
    
    @classmethod
    def create_new(cls, project_name: str = "Untitled Project") -> 'WorkflowSession':
        """Create a new workflow session"""
        timestamp = datetime.now().isoformat()
        session_id = f"session_{int(time.time())}"
        
        return cls(
            session_id=session_id,
            project_name=project_name,
            created_at=timestamp,
            last_updated=timestamp
        )


class WorkflowValidator:
    """Validates workflow transitions and data integrity"""
    
    @staticmethod
    def validate_sequence(sequence: str) -> tuple[bool, str]:
        """Validate protein sequence"""
        if not sequence or len(sequence.strip()) == 0:
            return False, "Sequence cannot be empty"
        
        clean_seq = sequence.strip().upper().replace(" ", "")
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = set(clean_seq) - valid_aa
        
        if invalid_chars:
            return False, f"Invalid amino acids: {', '.join(invalid_chars)}"
        
        if len(clean_seq) < 10:
            return False, "Sequence too short (minimum 10 residues)"
        
        if len(clean_seq) > 2000:
            return False, "Sequence too long (maximum 2000 residues)"
        
        return True, clean_seq
    
    @staticmethod
    def validate_pdb_content(pdb_content: str) -> tuple[bool, str]:
        """Validate PDB file content"""
        if not pdb_content:
            return False, "PDB content is empty"
        
        lines = pdb_content.split('\n')
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        
        if len(atom_lines) == 0:
            return False, "No ATOM records found in PDB"
        
        if len(atom_lines) < 10:
            return False, "Too few atoms in structure"
        
        return True, f"Valid PDB with {len(atom_lines)} atoms"
    
    @staticmethod
    def validate_binding_site_residues(residues: List[int], max_residue: int) -> tuple[bool, str]:
        """Validate binding site residue selection"""
        if not residues:
            return True, "No binding site specified (will analyze full interface)"
        
        if any(r < 1 or r > max_residue for r in residues):
            return False, f"Residue numbers must be between 1 and {max_residue}"
        
        return True, f"Valid binding site with {len(residues)} residues"
