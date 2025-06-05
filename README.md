# PlasAnn - Plasmid Annotation Tool

This tool is designed for annotating plasmid sequences from FASTA or GenBank files. It provides comprehensive annotation including coding sequences (CDS), origins of replication (oriC/oriV), origins of transfer (oriT), transposons, replicons, and non-coding RNAs.

## Features

- **Automatic Database Download**: Databases are downloaded automatically on first run
- **Flexible Input**: Single files or entire folders of FASTA/GenBank files
- **Gene Prediction**: Uses Prodigal for accurate gene calling
- **Mobile Element Detection**: Identifies transposons, insertion sequences, and origins
- **Visualizations**: Generates circular plasmid maps
- **Batch Processing**: Process multiple files simultaneously
- **Enhanced Annotation**: Optional UniProt BLAST for comprehensive protein annotation

## Installation

### Dependencies

This tool requires the following external programs:
- **BLAST+** (makeblastdb, blastn, blastx, blastp)
- **Prodigal** (CDS prediction)
- **Infernal** (cmscan, cmpress - for ncRNA detection)

### Install PlasAnn

```bash
pip install plasann
```

### Installing External Dependencies

**For macOS (Apple Silicon M1/M2):**
```bash
brew tap brewsci/bio
brew install blast prodigal infernal
```

**For Linux (Ubuntu/Debian):**
```bash
sudo apt install ncbi-blast+ prodigal infernal
```

**For other systems:**
Install from source following the official documentation for [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/), [Prodigal](https://github.com/hyattpd/Prodigal), and [Infernal](http://eddylab.org/infernal/).

### Verify Installation

```bash
PlasAnn --check-deps
```

## Usage

### Basic Commands

```bash
PlasAnn -i <input_file_or_directory> -o <output_directory> -t <file_type>
```

### Parameters

- `-i`, `--input`: Path to input file or directory containing FASTA or GenBank files
- `-o`, `--output`: Path to output directory where results will be stored
- `-t`, `--type`: Type of input files: `fasta`, `genbank`, or `auto`

### Examples

**Single FASTA file:**
```bash
PlasAnn -i plasmid.fasta -o results -t fasta
```

**Folder of FASTA files:**
```bash
PlasAnn -i fasta_folder/ -o results -t fasta
```

**GenBank file (retain existing annotations):**
```bash
PlasAnn -i plasmid.gb -o results -t genbank --retain
```

**GenBank file (re-annotate with Prodigal):**
```bash
PlasAnn -i plasmid.gb -o results -t genbank --overwrite
```

**Auto-detect mixed file types:**
```bash
PlasAnn -i mixed_folder/ -o results -t auto
```

**Enhanced annotation with UniProt:**
```bash
PlasAnn -i plasmid.fasta -o results -t fasta --uniprot-blast
```

### GenBank Processing Options

When using GenBank files, you can choose:
- `--retain`: Use existing CDS annotations in the GenBank file (default)
- `--overwrite`: Ignore existing annotations and re-predict genes using Prodigal

### Additional Options

- `--uniprot-blast`: Run enhanced UniProt BLAST annotation (slower but more comprehensive)
- `--min-identity`: Minimum identity percentage for UniProt BLAST hits (default: 50%)
- `--version`: Show version information
- `--check-deps`: Check external dependency status

## Output Files

For each input file, PlasAnn generates:

- **CSV Table**: Detailed annotation table with gene information
- **GenBank File**: Annotated GenBank file with all features
- **Plasmid Map**: Circular visualization of the annotated plasmid (PNG format)
- **Enhanced CSV**: Additional file with UniProt annotations (if `--uniprot-blast` used)

## Troubleshooting

**Check dependencies:**
```bash
PlasAnn --check-deps
```

**Common issues:**
- Ensure BLAST+, Prodigal, and Infernal are properly installed
- Check that input files are valid FASTA or GenBank format
- For very short sequences (<100bp), CDS prediction may not work well

## Contact

For questions or issues, contact: [hislam2@ur.rochester.edu](mailto:hislam2@ur.rochester.edu)

## Citation

If you use PlasAnn in your research, please cite:
```
PlasAnn: Comprehensive Plasmid Annotation Pipeline
Habibul Islam
University of Rochester
```
