# Plasmid Annotation Tool

This tool is designed for annotating large plasmid sequences from FASTA or GenBank files. It utilizes various databases and tools to provide a comprehensive annotation of plasmid sequences, including the detection of coding sequences (CDS), origins of replication, transposons, and more.

## Features

- When the program is run for the first time the Databases will get downloaded automatically. 
- The pipeline can take single fasta file or a folder of fasta files as input. 
- It can also take a genbank file or a folder of genbank files as input 


## Installation

### Dependencies

This script requires the following Python packages:
- `gdown`
- `argparse`
- `pandas`
- `biopython`
- `matplotlib`
- `pycirclize`

The pipelines can be installed using pip: 

```bash
pip install plasann
```

or conda:
```bash
conda install bioconda::PlasAnn
```

In arm64 osx Blast is not installable through conda. We suggest installing blast using [homebrew](https://formulae.brew.sh/formula/blast) \
Or you can also install command line blast from [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
## Usage

To run the pipeline, use the following command:


```bash
PlasAnn -i <input_file_or_directory> -o <output_directory> -t <file_type>
```

### Parameters

- `-i`, `--input`: Path to the input file or directory containing FASTA or GenBank files.
- `-o`, `--output`: Path to the output directory where results will be stored.
- `-t`, `--type`: Type of the input files, either `fasta` or `genbank`.

Upon choosing GenBank as the file type, you will be prompted to select one of the following options:
1. Annotate the existing CDS in the genbank file. 
2. Overwrite existing CDS in GenBank files using Prodigal.

### Outputs

- **CSV Tables:** Contain the annotation details for each plasmid.
- **GenBank Files:** Annotated GenBank files with updated feature annotations.
- **Plasmid Maps:** PNG images representing the annotated plasmid.

The scripts are uploaded in the Scripts folder in this repository

Regarding any issues [Contact me](hislam2@ur.rochester.edu)
