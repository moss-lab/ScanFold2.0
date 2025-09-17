# ScanFold2.0

ScanFold2.0 identifies RNA secondary structures with unusually stable folding by comparing native RNA sequences against randomized controls using machine learning-based z-score calculations.

## Installation

### Prerequisites
- Python 3.10 or higher
- Conda (recommended) or pip

### Quick Setup with Conda

```bash
# Clone repository
git clone https://github.com/moss-lab/ScanFold2.0.git
cd ScanFold2.0

# Create environment
conda env create -f environment.yml
conda activate scanfold2
```

### Alternative Setup with pip

#### macOS (Apple Silicon):
```bash
python3 -m venv scanfold_env
source scanfold_env/bin/activate
pip install -r requirements-mac.txt
```

#### Linux/Windows:
```bash
python3 -m venv scanfold_env
source scanfold_env/bin/activate  # Windows: scanfold_env\Scripts\activate
pip install -r requirements-linux.txt
```

**Note:** ViennaRNA must be installed separately via conda or system package manager.

## Input Requirements

- **Format**: FASTA file containing RNA sequence(s)
- **Sequence type**: RNA sequences (U's) or DNA (T's will be converted to U's)
- **Length**: Sequences should be longer than window size (default 120 nt)
- **Multiple sequences**: Each sequence in the FASTA will be analyzed separately

## Quick Start

```bash
# Basic usage
python ScanFold2.0.py input.fasta

# With custom parameters
python ScanFold2.0.py input.fasta -w 120 -s 1 -t 37 --shuffle di
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|----------|
| `-w` | Window size (nt) | 120 |
| `-s` | Step size (nt) | 1 |
| `-t` | Temperature (°C) | 37 |
| `--shuffle` | Shuffling type: `mono` or `di` | mono |
| `--folder` | Output folder name | ScanFold_YYYYMMDD_HHMMSS |

### Shuffling Types
- **`mono`**: Mononucleotide shuffling - preserves single nucleotide composition
- **`di`**: Dinucleotide shuffling - preserves dinucleotide frequencies for better background modeling

## Output Files

All results are saved to a timestamped folder (or custom named folder with `--folder`). Key outputs include:

- **`.tsv`**: Main results table with z-scores for each window
- **`.ct` files**: RNA structures in CT format (filtered at different z-score thresholds)
- **`.dbn` files**: Dot-bracket notation structures
- **`.wig` files**: Track files for genome browser visualization
- **`RESULTS_README.md`**: Detailed guide for viewing and interpreting results
- **`ScanFold_run.log`**: Complete analysis log

## Examples

### Basic RNA structure scanning
```bash
python ScanFold2.0.py my_rna.fasta
```

### High-resolution scanning with dinucleotide shuffling
```bash
python ScanFold2.0.py my_rna.fasta -w 80 -s 1 --shuffle di --folder my_analysis
```

### Temperature-specific folding
```bash
python ScanFold2.0.py viral_rna.fasta -t 30 -w 150
```

## Interpreting Results

- **Z-score < -1**: Potentially stable structures
- **Z-score < -2**: Highly stable, likely functional structures
- **Lower z-scores** indicate more stable folding compared to random sequences

## Visualizing in IGV

The WIG and BP files can be loaded into IGV (Integrative Genomics Viewer) for visualization:

1. **Download IGV** from https://igv.org if not already installed
2. **Load your reference genome** or sequence file
3. **Load ScanFold tracks**:
   - File → Load from File
   - Navigate to your ScanFold output folder
   - Select the WIG files:
     - `.scan-zscores.wig` - Z-score track
     - `.scan-MFE.wig` - Minimum free energy values
     - `.scan-ED.wig` - Ensemble diversity values
     - `.Zavg.wig` - Average Z-scores per base pair

4. **Customize visualization**:
   - Right-click tracks to adjust color schemes
   - Set z-score track to use blue-white-red gradient (blue for negative/stable)
   - Adjust track height for better visibility

## Troubleshooting

### Common Issues

**TensorFlow/Metal warnings on Mac**
- These are normal and can be ignored
- The tool will still run correctly

**Memory issues with large sequences**
- Try increasing step size (`-s`) to reduce number of windows
- Consider splitting very long sequences

**ViennaRNA not found**
- Ensure ViennaRNA is installed: `conda install -c bioconda viennarna`
- Check it's accessible: `RNAfold -h`

**NumPy version conflicts**
- Must use NumPy < 2.0 for TensorFlow compatibility
- Check with: `pip show numpy`

## Performance

- **Runtime**: ~1-2 seconds per 1000 nt (depending on step size)
- **Memory**: Generally < 2GB for typical sequences
- **GPU**: Will use Apple Metal on Mac or CUDA on Linux if available

## Citation

If you use ScanFold2.0, please cite:

Andrews RJ, Rouse WB, O'Leary CA, Booher NJ, Moss WN. (2023). ScanFold 2.0: a rapid approach for identifying potential structured RNA targets in genomes and transcriptomes. *PeerJ* 11:e16061. https://doi.org/10.7717/peerj.16061

## Support

For issues or questions, please open an issue on GitHub