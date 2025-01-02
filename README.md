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

The pipelines can be installed using pip using 

```bash
pip install plasann
```

and using conda:
```bash
conda install bioconda::PlasAnn
```

### Installing Prodigal

You can download and install Prodigal from its official repository: from [here](https://github.com/hyattpd/Prodigal/wiki/installation) \


### Installing command line blast

Command line blast needs to be installed. Follow the installation instruction [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/) \
Remember to append the path to the new BLAST bin directory to the existing PATH setting. \

For linux, Windows and osx64 users Blast can be installed using : conda install bioconda::blast\
For osxARM64 we suggest we suggest usiing Blast from the link above Due to compatibility issues, . \

## Usage

To run the pipeline, use the following command:


```bash
PlasAnn -i <input_file_or_directory> -o <output_directory> -t <file_type>
```

### Parameters

- `-i`, `--input`: Path to the input file or directory containing FASTA or GenBank files.
- `-o`, `--output`: Path to the output directory where results will be stored.
- `-t`, `--type`: Type of the input files, either `fasta` or `genbank`.

### Example Commands

#### For FASTA Files:
```bash
python annotate_plasmids.py -i /path/to/fasta/files -o /path/to/output/directory -t fasta
```

#### For GenBank Files:
```bash
python annotate_plasmids.py -i /path/to/genbank/files -o /path/to/output/directory -t genbank
```

Upon choosing GenBank as the file type, you will be prompted to select one of the following options:
1. Annotate the existing CDS in the genbank file. 
2. Overwrite existing CDS in GenBank files using Prodigal.

### Outputs

- **CSV Tables:** Contain the annotation details for each plasmid.
- **GenBank Files:** Annotated GenBank files with updated feature annotations.
- **Plasmid Maps:** PNG images representing the annotated plasmid.

The scripts are uploaded in the Scripts folder in this repository

Regarding any issues [Contact me](hislam2@ur.rochester.edu)