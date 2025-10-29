# Protein Structure Prediction with NVIDIA AI

A Streamlit web application that predicts protein 3D structures from amino acid sequences using NVIDIA's Cloud Functions with multiple AI models including AlphaFold2, OpenFold2, and Boltz2.

## 🌟 Features

- 🧬 **Multiple AI Models**: Choose from OpenFold2, AlphaFold2, AlphaFold2 Multimer, and Boltz2
- 🎮 **Interactive 3D Visualization**: View predicted structures with py3Dmol
- 📥 **Download PDB Files**: Save predicted structures as standard PDB files
- ✅ **Sequence Validation**: Automatic validation of amino acid sequences
- 🎯 **Demo Mode**: Test the interface with mock predictions
- 📊 **Robust API Handling**: Multiple payload formats and error handling
- 🔄 **Asynchronous Processing**: Handles long-running predictions with progress tracking

## 🚀 Quick Start

### Option 1: Use the launcher script
```bash
./launch.sh
```

### Option 2: Manual launch
```bash
# Activate virtual environment
source .venv/bin/activate

# Install dependencies (if needed)
pip install -r requirements.txt

# Run the application
streamlit run app_v2.py --server.port 8502
```

### Option 3: Run both versions
- **Original version** (port 8501): `streamlit run app.py`
- **Improved version** (port 8502): `streamlit run app_v2.py --server.port 8502`

## 🌐 Access the Application

Once running, open your browser and navigate to:
- **Main App**: http://localhost:8502
- **Original App**: http://localhost:8501 (if running)

## 🧬 Available Models

| Model | Description | Best For |
|-------|-------------|----------|
| **OpenFold2** | Open source protein structure prediction | General protein structures |
| **AlphaFold2** | DeepMind's protein structure prediction | High accuracy single chains |
| **AlphaFold2 Multimer** | AlphaFold2 for protein complexes | Multi-chain proteins |
| **Boltz2** | Advanced protein structure prediction | Latest improvements |

## 📝 Usage Guide

1. **Select a Model**: Choose from the available AI models in the sidebar
2. **Enter API Key**: Input your NVIDIA API key (pre-configured)
3. **Input Sequence**: 
   - Paste an amino acid sequence (10-2000 residues)
   - Or select from provided examples
   - Use standard single-letter amino acid codes
4. **Predict Structure**: Click the "Predict Structure" button
5. **View Results**: 
   - Interactive 3D visualization
   - Download PDB file
   - View raw PDB content

## 🧪 Demo Mode

Enable "Demo Mode" in the sidebar to test the interface without making real API calls. This generates mock PDB structures for demonstration purposes.

## 📋 Example Sequences

The application includes several pre-loaded examples:
- **Insulin B-chain** (30 AA): Short peptide hormone
- **Lysozyme fragment** (140 AA): Antimicrobial enzyme
- **Cytochrome C fragment** (50 AA): Electron transport protein
- **Sample sequence** (27 AA): Simple test sequence

## 🔧 Configuration

### Environment Variables
```bash
export NVIDIA_API_KEY="your-api-key-here"
```

### API Configuration
- **Endpoint**: `https://api.nvcf.nvidia.com/v2/nvcf/pexec/functions/{model_id}`
- **Authentication**: Bearer token with NVIDIA API key
- **Timeout**: 300 seconds for structure prediction
- **Polling**: Automatic handling of asynchronous requests

## 🛠 Technical Details

### Architecture
- **Frontend**: Streamlit for web interface
- **API Client**: Direct REST API calls to NVIDIA Cloud Functions
- **Visualization**: py3Dmol for 3D molecular rendering
- **File Handling**: Standard PDB format output
- **Error Handling**: Robust error handling and fallback mechanisms

### File Structure
```
nvidia/
├── app.py                 # Original application
├── app_v2.py             # Improved application (recommended)
├── protein_models.py     # Model configuration
├── requirements.txt      # Python dependencies
├── launch.sh            # Launch script
├── setup.sh             # Setup script
├── test_app.py          # Test suite
├── README.md            # This file
├── QUICK_START.md       # Quick start guide
└── .env.template        # Environment variables template
```

### Key Improvements in v2
- Multiple model support with easy selection
- Robust API error handling with multiple payload formats
- Better progress tracking for long-running predictions
- Demo mode for testing without API calls
- Enhanced user interface with metrics and better organization
- Comprehensive error messages and troubleshooting

## 🔍 Troubleshooting

### Common Issues

1. **API Key Issues**
   - Ensure your NVIDIA API key is valid and active
   - Check that the key has access to Cloud Functions
   - Verify the key in the sidebar configuration

2. **Model Selection**
   - Different models may have different input requirements
   - Try switching models if one fails
   - Check model status in the NVIDIA console

3. **Sequence Validation**
   - Use only standard amino acid codes (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
   - Ensure sequences are 10-2000 amino acids long
   - Remove any special characters or spaces

4. **Long Processing Times**
   - Structure prediction can take several minutes
   - The app shows progress for asynchronous requests
   - Check your internet connection if requests timeout

5. **Visualization Issues**
   - If 3D visualization doesn't load, check browser JavaScript settings
   - Try downloading the PDB file and using external visualization tools
   - Check the raw PDB content in the expandable section

### Error Messages

| Error | Solution |
|-------|----------|
| "Invalid amino acid characters" | Remove non-standard characters from sequence |
| "Sequence too short/long" | Adjust sequence length to 10-2000 amino acids |
| "API Error: 401" | Check API key validity and permissions |
| "API Error: 429" | Rate limit exceeded, wait and retry |
| "All payload formats failed" | Try a different model or enable demo mode |

## 📊 Performance

- **Sequence Length**: Optimized for sequences up to 2000 amino acids
- **Processing Time**: 2-10 minutes depending on sequence length and model
- **Models**: All models run on NVIDIA's cloud infrastructure
- **Rate Limits**: Subject to NVIDIA Cloud Functions limits

## 🔐 Security

- API keys are handled securely and not logged
- All communication uses HTTPS
- Sequences are processed by NVIDIA's secure cloud infrastructure
- No data is stored permanently on our servers

## 🆘 Support

For issues related to:
- **NVIDIA Cloud Functions**: Contact NVIDIA support or check their status page
- **This Application**: Check the troubleshooting section or GitHub issues
- **Streamlit**: Refer to Streamlit documentation
- **Protein Biology**: Consult relevant scientific literature

## 📄 License

This project is for educational and research purposes. Please check NVIDIA's terms of service for commercial usage of their Cloud Functions.

## 🙏 Acknowledgments

- NVIDIA for providing Cloud Functions and AI models
- DeepMind for AlphaFold2
- OpenFold team for OpenFold2
- py3Dmol for molecular visualization
- Streamlit for the web framework
