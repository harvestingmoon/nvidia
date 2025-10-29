# Binding Protein Design Workflow

A comprehensive web application for designing and analyzing protein binders using NVIDIA's AI models and a multi-step workflow: **ESMFold → RFDiffusion → ProteinMPNN → DiffDock**

## 🎯 Features

### Multi-Step Workflow
- **Stage 1: Target Input** - Define your target protein (sequence, PDB file, or PDB ID)
- **Stage 2: Target Prediction** - Predict 3D structure using ESMFold/AlphaFold2
- **Stage 3: Binder Design** - Design or input your binding protein sequence
- **Stage 4: Binder Prediction** - Predict binder 3D structure
- **Stage 5: Complex Analysis** - Analyze binding interface and quality
- **Stage 6: Results** - View comprehensive results and export data

### State Management
- **Session Persistence** - Save and load entire workflow sessions
- **Progress Tracking** - Visual progress stepper showing completed/pending stages
- **Data Validation** - Automatic validation at each step
- **Checkpoint System** - Resume from any stage
- **Project Notes** - Add notes and documentation

### Binding Analysis
- **Interface Detection** - Identify residues at binding interface
- **Quality Scoring** - Automated quality assessment (0-100 score)
- **Contact Analysis** - Detailed contact pairs and distances
- **3D Visualization** - Interactive structure viewer with Mol*
- **Export Options** - Download PDB files and analysis reports

## 🚀 Quick Start

### Option 1: Run the Binding Workflow App

```bash
# Activate virtual environment
source .venv/bin/activate  # macOS/Linux
# or
.venv\Scripts\activate  # Windows

# Run the binding workflow app
streamlit run binding_workflow_app.py
```

### Option 2: Run the Original Single-Structure App

```bash
streamlit run app_v2.py
```

## 📋 Prerequisites

1. **Python 3.10+**
2. **NVIDIA API Key** (optional - demo mode available)
   - Get your API key from [NVIDIA NGC](https://catalog.ngc.nvidia.com/)
3. **Dependencies** (install via requirements.txt)

```bash
pip install -r requirements.txt
```

## 🔧 Configuration

### API Setup
1. Enable demo mode in the sidebar, OR
2. Enter your NVIDIA API key in the sidebar

### Workflow Configuration
- **Interface Cutoff**: Adjust distance threshold for interface detection (default: 5.0 Å)
- **Binding Site**: Optionally specify target residues for focused analysis
- **Model Selection**: Choose between OpenFold2 (faster) or AlphaFold2 (more accurate)

## 📊 Workflow Details

### Target Input Stage
- **Input Types**:
  - Amino acid sequence (10-2000 residues)
  - PDB file upload
  - PDB ID from RCSB PDB
- **Validation**: Automatic sequence validation
- **Binding Site Specification**: Optional residue selection

### Structure Prediction
- **Models Available**:
  - OpenFold2 (fast, ~2-5 minutes)
  - AlphaFold2 (accurate, ~5-10 minutes)
  - AlphaFold2 Multimer
  - Boltz2
- **Progress Tracking**: Real-time status updates
- **Retry Logic**: Automatic retry with multiple payload formats

### Binder Design
- **Manual Sequence**: Enter designed sequence directly
- **RFDiffusion**: Coming soon - generative backbone design
- **Template-based**: Coming soon - design from existing scaffolds
- **Sequence Analysis**: Hydrophobic, charged residue composition

### Complex Analysis
- **Docking Methods**:
  - Simple overlay (immediate)
  - DiffDock (coming soon - ML-based docking)
- **Interface Analysis**:
  - Residue contact detection
  - Distance calculations
  - Buried surface area estimation
- **Quality Metrics**:
  - Number of contacts
  - Average/minimum distances
  - Steric clash detection
  - Overall quality score (0-100)

## 📁 Project Structure

```
nvidia/
├── binding_workflow_app.py      # Main workflow application
├── workflow_state.py             # State management and data structures
├── binding_analysis.py           # Interface analysis algorithms
├── app_v2.py                     # Original single-structure app
├── protein_models.py             # Model configurations
├── pdb_viewer.py                 # 3D visualization utilities
├── requirements.txt              # Python dependencies
└── README_WORKFLOW.md            # This file
```

## 🔬 Analysis Outputs

### Quality Scoring
- **Grade A (85-100)**: Excellent binding interface, ready for validation
- **Grade B (70-84)**: Good interface, minor optimization suggested
- **Grade C (50-69)**: Moderate quality, optimization recommended
- **Grade D (30-49)**: Poor interface, redesign needed
- **Grade F (<30)**: Insufficient interface, complete redesign required

### Scoring Criteria
1. **Number of Contacts** (40 points)
   - ≥15 contacts: 40 pts
   - 10-14 contacts: 30 pts
   - 5-9 contacts: 15 pts
   - <5 contacts: 0 pts

2. **Average Distance** (40 points)
   - 3.5-4.5 Å: 40 pts
   - 3.0-5.0 Å: 25 pts
   - <3.0 Å or >5.0 Å: 10-15 pts

3. **Steric Clashes** (20 points)
   - Min distance ≥2.8 Å: 20 pts
   - Min distance ≥2.5 Å: 10 pts
   - Min distance <2.5 Å: 0 pts

## 💾 Session Management

### Save Session
```json
{
  "session_id": "session_1698765432",
  "project_name": "MyBinderProject",
  "created_at": "2025-10-30T10:30:00",
  "current_stage": "complex_analysis",
  "target": {...},
  "binder": {...},
  "complex": {...}
}
```

### Load Session
- Upload saved JSON file
- Resume from any stage
- All data and progress preserved

## 🎨 UI Features

### Progress Visualization
- Visual stepper showing all 6 stages
- Color-coded status (complete, in-progress, pending, failed)
- Progress bar showing overall completion

### Interactive Elements
- 3D structure viewer with mouse controls
- Collapsible advanced options
- Real-time validation feedback
- Downloadable results at each stage

### Responsive Design
- Wide layout for complex visualizations
- Column-based metric displays
- Expandable sections for detailed info

## 📤 Export Options

### Individual Files
- Target PDB
- Binder PDB
- Complex PDB
- Session JSON

### Complete Project Export
- All structures
- Analysis results
- Workflow metadata
- Project notes

## 🔮 Coming Soon

### RFDiffusion Integration
- Generative backbone design
- Motif scaffolding
- Symmetric designs

### ProteinMPNN Integration
- Sequence design for backbones
- Multiple sequence suggestions
- Sequence scoring

### DiffDock Integration
- ML-based molecular docking
- Multiple binding poses
- Confidence scoring
- Small molecule ligand support

### Enhanced Analysis
- Binding affinity prediction
- Rosetta energy calculations
- Molecular dynamics setup
- Mutation suggestions

## 🐛 Troubleshooting

### API Issues
- **Timeout**: Try demo mode or shorter sequences
- **Auth Error**: Verify API key is correct
- **Rate Limit**: Wait a few minutes and retry

### Structure Issues
- **No ATOM records**: Check PDB file format
- **Visualization error**: Try downloading and viewing in PyMOL
- **Missing CA atoms**: Structure may be incomplete

### Interface Issues
- **No contacts found**: Adjust distance cutoff or check alignment
- **Low quality score**: Review sequence design and interface residues

## 📚 Resources

- [NVIDIA NGC Catalog](https://catalog.ngc.nvidia.com/)
- [AlphaFold2 Paper](https://www.nature.com/articles/s41586-021-03819-2)
- [RFDiffusion Paper](https://www.nature.com/articles/s41586-023-06415-8)
- [ProteinMPNN Paper](https://www.science.org/doi/10.1126/science.add2187)
- [DiffDock Paper](https://arxiv.org/abs/2210.01776)

## 📝 Citation

If you use this workflow in your research, please cite:
- NVIDIA NIM APIs
- AlphaFold2 (Jumper et al., 2021)
- OpenFold (Ahdritz et al., 2022)

## 🤝 Contributing

This is a demonstration application. For production use:
1. Add proper error handling
2. Implement authentication
3. Add result validation
4. Include experimental verification workflows

## 📄 License

This project uses NVIDIA's APIs and is subject to their terms of service.

## 💡 Tips

1. **Start with Demo Mode** to understand the workflow
2. **Use OpenFold2** for faster predictions during testing
3. **Save sessions frequently** to preserve progress
4. **Export results** before starting new designs
5. **Review interface residues** for optimization targets

## 🎓 Learning Resources

### Protein Structure Prediction
- Understanding pLDDT confidence scores
- Interpreting predicted structures
- Validation methods

### Binding Interface Design
- Principles of protein-protein interactions
- Hotspot residues
- Interface optimization strategies

### Computational Drug Design
- Structure-based design
- Virtual screening
- Experimental validation

---

**Ready to design your first binder?** Run `streamlit run binding_workflow_app.py` and start exploring! 🚀
