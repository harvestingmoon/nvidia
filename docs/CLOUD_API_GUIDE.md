# NVIDIA AlphaFold2 Multimer - Cloud API

This tool predicts protein structures using NVIDIA's hosted AlphaFold2 and AlphaFold2-Multimer APIs.

## 🚀 Quick Start

### Prerequisites
```bash
# Ensure NGC API key is set
export NGC_CLI_API_KEY='your_key_here'

# Get your key from: https://org.ngc.nvidia.com/setup/api-key
```

### Installation
```bash
# Install required packages
pip install requests
```

## 📖 Usage

### Single Protein (Monomer)

Predict the structure of a single protein chain:

```bash
# From sequence
.venv/bin/python nvidia_cloud_api.py "MKFLKFSLLTAVLLSVVFAFSSCGDDDDTYPYDVPDYA" -o protein.pdb

# From FASTA file
.venv/bin/python nvidia_cloud_api.py examples/monomer_example.fasta -o protein.pdb
```

### Protein Complex (Multimer)

Predict the structure of a protein complex with multiple chains:

```bash
# Two chains - direct input
.venv/bin/python nvidia_cloud_api.py \
    "MKFLKFSLLTAVLLSVVFA" \
    "GAGYPYDVPDYAGGGEQKL" \
    -o complex.pdb

# Multiple chains from FASTA file
.venv/bin/python nvidia_cloud_api.py examples/multimer_example.fasta -o complex.pdb

# Three or more chains
.venv/bin/python nvidia_cloud_api.py \
    "SEQUENCE1" \
    "SEQUENCE2" \
    "SEQUENCE3" \
    -o three_chain_complex.pdb
```

### FASTA Format for Multimer

Create a FASTA file with multiple sequences:

```fasta
>Chain_A
MKFLKFSLLTAVLLSVVFAFSSCGDDDDTYPYDVPDYA
>Chain_B
GAGYPYDVPDYAGGGEQKLISEEDLLRKRREQLKHKLE
>Chain_C
AEAEDLLRKRREQLKHKLEQGGFGFGFGQWDWERVMG
```

## 🎯 Command Options

```bash
usage: nvidia_cloud_api.py [-h] [-o OUTPUT] [-a {mmseqs2,jackhmmer}] 
                           [-t TIMEOUT] [--monomer] [--json]
                           sequences [sequences ...]

Arguments:
  sequences              One or more sequences or FASTA file path
  
Options:
  -o, --output FILE      Output PDB filename (default: predicted_structure.pdb)
  -a, --algorithm NAME   MSA algorithm: mmseqs2 or jackhmmer (default: mmseqs2)
  -t, --timeout SECONDS  Request timeout (default: 600)
  --monomer              Force monomer prediction (use only first sequence)
  --json                 Also save full JSON response with confidence scores
  -h, --help             Show help message
```

## 📋 Examples

### Example 1: Simple Monomer
```bash
.venv/bin/python nvidia_cloud_api.py examples/monomer_example.fasta -o outputs/protein.pdb
```

### Example 2: Antibody-Antigen Complex
```bash
# Heavy chain + Light chain + Antigen
.venv/bin/python nvidia_cloud_api.py \
    examples/heavy_chain.fasta \
    examples/light_chain.fasta \
    examples/antigen.fasta \
    -o outputs/antibody_complex.pdb
```

### Example 3: Homodimer
```bash
# Same sequence twice for homodimer
.venv/bin/python nvidia_cloud_api.py \
    "MKFLKFSLLTAVLLSVVFA" \
    "MKFLKFSLLTAVLLSVVFA" \
    -o outputs/homodimer.pdb
```

### Example 4: With Confidence Data
```bash
# Save JSON with pLDDT and PAE scores
.venv/bin/python nvidia_cloud_api.py \
    examples/multimer_example.fasta \
    -o outputs/complex.pdb \
    --json
```

## 🔍 Understanding Output

### PDB File
The `.pdb` file contains the predicted 3D structure:
- For monomers: Single chain labeled 'A'
- For multimers: Multiple chains labeled 'A', 'B', 'C', etc.
- Can be opened in PyMOL, ChimeraX, or Mol*

### JSON File (with --json flag)
Contains additional data:
- `pdb`: The structure in PDB format
- `plddt`: Per-residue confidence scores (0-100)
- `pae`: Predicted aligned error matrix (for multimers)
- `ranking_confidence`: Overall model quality score
- `_metadata`: Chain information

## 🧬 Multimer vs Monomer

| Feature | Monomer | Multimer |
|---------|---------|----------|
| Input | Single sequence | Multiple sequences |
| Use case | Single protein | Protein complexes |
| Endpoint | `/alphafold2` | `/alphafold2-multimer` |
| Output | Single chain PDB | Multi-chain PDB |
| Interface prediction | ❌ | ✅ |
| PAE scores | ❌ | ✅ |

## 🎨 Visualization

### Using PyMOL
```bash
pymol outputs/complex.pdb
```

### Using ChimeraX
```bash
chimerax outputs/complex.pdb
```

### Using Streamlit App
```bash
bash scripts/launch.sh
```

## ⚠️ Important Notes

1. **Sequence Length Limits**:
   - Monomer: Up to ~2,000 amino acids
   - Multimer: Total length typically limited to ~3,000 amino acids
   - Longer sequences may timeout or fail

2. **Multimer Format**:
   - Sequences are joined with `:` separator internally
   - Order matters for interface prediction
   - Each sequence must be valid protein sequence (20 amino acids)

3. **Computation Time**:
   - Monomer: 30 seconds - 5 minutes
   - Multimer: 2 - 20 minutes depending on size
   - Large complexes may require longer timeout: `-t 1200`

4. **Valid Amino Acids**:
   ```
   A C D E F G H I K L M N P Q R S T V W Y
   ```
   No special characters or lowercase allowed.

## 🔧 Troubleshooting

### "NGC_CLI_API_KEY not found"
```bash
export NGC_CLI_API_KEY='nvapi-YOUR_KEY_HERE'
```

### "Request timed out"
Increase timeout for large complexes:
```bash
.venv/bin/python nvidia_cloud_api.py sequences.fasta -o output.pdb -t 1200
```

### "Invalid amino acids"
Check your sequence only contains valid amino acids (A-Y, no B/J/O/U/X/Z)

### API Rate Limits
If you hit rate limits, wait a few minutes and try again.

## 📚 Additional Resources

- [NVIDIA NIM Documentation](https://docs.nvidia.com/nim/)
- [AlphaFold2 Paper](https://www.nature.com/articles/s41586-021-03819-2)
- [AlphaFold-Multimer Paper](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v2)

## 🎯 Use Cases

### Drug Discovery
- Predict antibody-antigen binding
- Model protein-protein interactions
- Screen potential drug targets

### Structural Biology
- Predict protein complex structures
- Study protein-protein interfaces
- Validate experimental structures

### Protein Engineering
- Design protein complexes
- Optimize binding interfaces
- Engineer stable dimers/trimers

## 🚀 Run All Examples

```bash
bash examples/run_examples.sh
```

This will run through all example use cases interactively.
