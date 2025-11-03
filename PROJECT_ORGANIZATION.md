
nvidia/
│
├── frontend/                      # Web Applications
│   ├── app_v2.py                  # Single structure prediction app
│   ├── binding_workflow_app.py    # Multi-step binding workflow app
│   └── __init__.py                # Package initialization
│
├── workflow/                      # Workflow Management
│   ├── workflow_state.py          # State management & sessions
│   ├── binding_analysis.py        # Interface analysis algorithms
│   └── __init__.py                # Package exports
│
├── core/                          # Core Utilities
│   ├── protein_models.py          # Model configurations (NVIDIA APIs)
│   ├── pdb_viewer.py              # 3D visualization utilities
│   └── __init__.py                # Package exports
│
├── scripts/                       # Launch Scripts
│   ├── launch.sh                  # Launch single structure app (port 8502)
│   ├── launch_workflow.sh         # Launch workflow app (port 8501)
│   └── setup.sh                   # Environment setup script
│
├── docs/                          # Documentation
│   ├── README.md                  # Main documentation
│   ├── README_WORKFLOW.md         # Workflow app guide
│   ├── APP_COMPARISON.md          # App comparison
│   ├── QUICK_START.md             # Quick start guide
│   ├── ALPHAFOLD2_IMPROVEMENTS.md # AlphaFold2 notes
│   └── PDB_VIEWING_IMPROVEMENTS.md# PDB viewer notes
│
├── Configuration Files
│   ├── requirements.txt           # Python dependencies
│   ├── .env.template              # Environment template
│   └── README.md                  # Main project README
│
└── Other
    ├── old_code/                  # Archived code
    ├── .venv/                     # Virtual environment
    ├── .git/                      # Git repository
    └── __pycache__/               # Python cache
```

## Import Changes

### Before (Old Structure)
```python
from protein_models import PROTEIN_MODELS
from workflow_state import WorkflowSession
from binding_analysis import parse_pdb_content
```

### After (New Structure)
```python
from core.protein_models import PROTEIN_MODELS
from workflow.workflow_state import WorkflowSession
from workflow.binding_analysis import parse_pdb_content
```

## Running Applications

### Method 1: Using Scripts (Recommended)

```bash
# From project root
cd /Users/bytedance/Documents/nvidia

# Single structure prediction app
bash scripts/launch.sh

# Binding workflow app
bash scripts/launch_workflow.sh
```

### Method 2: Direct Command

```bash
# Activate environment
source .venv/bin/activate

# Set PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Run single structure app
streamlit run frontend/app_v2.py --server.port 8502

# Run workflow app
streamlit run frontend/binding_workflow_app.py --server.port 8501
```

## Key Benefits

### `frontend/`
- **Purpose**: User-facing applications
- **Contains**: Streamlit web apps
- **Imports from**: `workflow/`, `core/`

### `workflow/`
- **Purpose**: Business logic and state management
- **Contains**: Workflow engine, analysis algorithms
- **Imports from**: `core/`
- **Imported by**: `frontend/`

### `core/`
- **Purpose**: Shared utilities
- **Contains**: Model configs, visualization tools
- **Imports from**: Standard libraries
- **Imported by**: `frontend/`, `workflow/`

### `scripts/`
- **Purpose**: Automation
- **Contains**: Launch and setup scripts
- **Sets**: Environment variables (PYTHONPATH)

### `docs/`
- **Purpose**: Documentation
- **Contains**: All markdown documentation files
- **No code**: Pure documentation

## Bug Fixes Applied

### 1. Input Type Mapping Issue
**Problem**: ValueError when loading saved sessions
```python
# Before (buggy)
target.input_type.replace("_", " ").title()  # "pdb_file" → "Pdb File" ❌

# After (fixed)
input_type_map = {
    "sequence": "Sequence",
    "pdb_file": "PDB File",  # Correct mapping ✅
    "pdb_id": "PDB ID"
}
```

### 2. Import Path Updates
**Problem**: Module not found errors
```python
# Before
from app_v2 import generate_mock_pdb  # ❌

# After
from frontend.app_v2 import generate_mock_pdb  # ✅
```

### 3. PYTHONPATH Configuration
**Problem**: Python couldn't find modules
```bash
# Added to launch scripts
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```


### To Run the Apps:
```bash
# Navigate to project root
cd /Users/bytedance/Documents/nvidia

# Run workflow app
bash scripts/launch_workflow.sh
```

### To Develop:
1. Edit files in their respective folders
2. Frontend changes → `frontend/`
3. Workflow logic → `workflow/`
4. Shared utilities → `core/`
5. Documentation → `docs/`

## Tips

1. **Always run from project root**: Scripts expect to be run from `/Users/bytedance/Documents/nvidia`
2. **Use the launch scripts**: They handle PYTHONPATH and environment setup
3. **Check import paths**: Use full module paths (e.g., `from workflow.xxx import yyy`)
4. **Save sessions**: Workflow app supports save/load - use it!

## Quick Reference

| Task | Command |
|------|---------|
| Launch workflow app | `bash scripts/launch_workflow.sh` |
| Launch single app | `bash scripts/launch.sh` |
| Setup environment | `bash scripts/setup.sh` |
| Install dependencies | `pip install -r requirements.txt` |
| View docs | `cat docs/README_WORKFLOW.md` |

---

**Everything is now organized and ready to use!**
