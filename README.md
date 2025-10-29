# Protein Binding Design Workflow

A comprehensive platform for designing and analyzing protein binders using NVIDIA's AI models.

## 🎯 Project Overview

This project provides two powerful applications:
1. **Single Structure Prediction** - Fast protein structure prediction
2. **Binding Workflow** - Complete binding protein design pipeline

**Workflow Pipeline**: ESMFold → RFDiffusion → ProteinMPNN → DiffDock

## 📁 Project Structure

```
nvidia/
├── frontend/                    # Web applications
│   ├── app_v2.py               # Single structure prediction app
│   ├── binding_workflow_app.py # Multi-step binding workflow app
│   └── __init__.py
│
├── workflow/                    # Workflow management
│   ├── workflow_state.py       # State management & data structures
│   ├── binding_analysis.py     # Interface analysis algorithms
│   └── __init__.py
│
├── core/                        # Core utilities
│   ├── protein_models.py       # Model configurations
│   ├── pdb_viewer.py          # 3D visualization utilities
│   └── __init__.py
│
├── scripts/                     # Launch scripts
│   ├── launch.sh              # Launch single structure app
│   ├── launch_workflow.sh     # Launch workflow app
│   └── setup.sh               # Environment setup
│
├── docs/                        # Documentation
│   ├── README_WORKFLOW.md     # Workflow app documentation
│   ├── APP_COMPARISON.md      # App comparison guide
│   ├── QUICK_START.md         # Quick start guide
│   └── *.md                   # Other documentation
│
├── old_code/                    # Archived code
├── requirements.txt             # Python dependencies
├── .env.template               # Environment template
└── README.md                   # This file
```

## 🚀 Quick Start

### 1. Setup Environment

```bash
# Clone or navigate to project
cd nvidia

# Run setup script
bash scripts/setup.sh
```

### 2. Run Applications

**Option A: Single Structure Prediction**
```bash
bash scripts/launch.sh
# Opens at http://localhost:8502
```

**Option B: Binding Workflow**
```bash
bash scripts/launch_workflow.sh
# Opens at http://localhost:8501
```

### 3. Configure API

- Get your NVIDIA API key from [NVIDIA NGC](https://catalog.ngc.nvidia.com/)
- Enter it in the app sidebar, OR
- Enable "Demo Mode" to test without API

## 📚 Documentation

- **[Workflow Guide](docs/README_WORKFLOW.md)** - Complete workflow documentation
- **[App Comparison](docs/APP_COMPARISON.md)** - Which app to use?
- **[Quick Start](docs/QUICK_START.md)** - Get started quickly

## 🔬 Features

### Single Structure App
- ✅ Fast structure prediction
- ✅ Multiple model support (ESMFold, AlphaFold2, Boltz2)
- ✅ 3D visualization
- ✅ Confidence scores (pLDDT)
- ✅ PDB export

### Binding Workflow App
- ✅ 6-stage workflow pipeline
- ✅ Target & binder structure prediction
- ✅ Interface analysis
- ✅ Quality scoring (0-100)
- ✅ Session save/load
- ✅ Project management
- ✅ Complex visualization
- ✅ Comprehensive export

## 🛠️ Technology Stack

**Backend:**
- Python 3.10+
- NVIDIA NIM APIs
- NumPy for calculations
- BioPython for structure handling

**Frontend:**
- Streamlit for web UI
- Plotly for interactive plots
- py3Dmol for 3D visualization

**Workflow:**
- State machine pattern
- Session persistence
- Automated validation

## 📊 Workflow Stages

1. **Target Input** - Define target protein
2. **Target Prediction** - Predict 3D structure
3. **Binder Design** - Design binder sequence
4. **Binder Prediction** - Predict binder structure
5. **Complex Analysis** - Analyze binding interface
6. **Results** - View and export complete analysis

## 🔧 Development

### Project Organization

- **`frontend/`** - User-facing Streamlit applications
- **`workflow/`** - State management and analysis logic
- **`core/`** - Shared utilities and configurations
- **`scripts/`** - Automation and launch scripts
- **`docs/`** - All documentation files

### Adding New Features

1. **Frontend changes** → Edit files in `frontend/`
2. **Workflow logic** → Edit files in `workflow/`
3. **Shared utilities** → Edit files in `core/`
4. **Documentation** → Add to `docs/`

### Import Structure

```python
# From frontend apps
from workflow.workflow_state import WorkflowSession
from workflow.binding_analysis import parse_pdb_content
from core.protein_models import PROTEIN_MODELS

# From workflow modules
from workflow import WorkflowSession, find_interface_residues
from core import PROTEIN_MODELS
```

## 📦 Dependencies

All dependencies are in `requirements.txt`:
- streamlit >= 1.28.0
- requests >= 2.31.0
- plotly >= 5.17.0
- numpy >= 1.24.0
- biopython >= 1.79
- pandas >= 1.5.0
- py3Dmol >= 2.0.0
- And more...

Install with:
```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## 🎓 Use Cases

### Research
- Design novel protein binders
- Analyze protein-protein interactions
- Validate computational designs

### Education
- Learn structure prediction
- Understand binding interfaces
- Practice computational biology

### Drug Discovery
- Design therapeutic proteins
- Screen binding candidates
- Optimize binding affinity

## 🤝 Contributing

This is a demonstration project. For production use:
1. Add comprehensive error handling
2. Implement authentication
3. Add result validation
4. Include experimental verification workflows

## 📄 License

This project uses NVIDIA's APIs and is subject to their terms of service.

## 🔗 Resources

- [NVIDIA NGC Catalog](https://catalog.ngc.nvidia.com/)
- [AlphaFold2 Paper](https://www.nature.com/articles/s41586-021-03819-2)
- [RFDiffusion Paper](https://www.nature.com/articles/s41586-023-06415-8)
- [Streamlit Documentation](https://docs.streamlit.io/)

## 💡 Tips

1. **Start with demo mode** to explore features
2. **Use single structure app** for quick predictions
3. **Use workflow app** for complete design projects
4. **Save sessions frequently** to preserve work
5. **Export results** before starting new projects

## 🐛 Troubleshooting

### Import Errors
If you see import errors, ensure you're running from the project root:
```bash
cd /path/to/nvidia
bash scripts/launch_workflow.sh
```

### API Issues
- Verify API key is correct
- Try demo mode to test interface
- Check NVIDIA service status

### Structure Issues
- Validate PDB format
- Check minimum sequence length (10 AA)
- Try different prediction models

## 📞 Support

- Check documentation in `docs/` folder
- Review code comments
- Test with demo mode first

---

**Ready to start?** Run `bash scripts/launch_workflow.sh` to begin! 🚀
