# PlasAnn
Plasmid annotation pipeline

This tool is designed for annotating plasmid sequences from FASTA or GenBank files. It utilizes various databases and tools to provide a comprehensive annotation of plasmid sequences, including the detection of coding sequences (CDS), origins of replication, transposons, and more.

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

Additionally, it relies on external tools such as Prodigal. Ensure these dependencies are installed:

```bash
pip install gdown pandas biopython matplotlib
```

### Installing Prodigal

You can download and install Prodigal from its official repository: from [here](https://github.com/hyattpd/Prodigal/wiki/installation)

### Installing command line blast

Command line blast needs to be installed. Follow the installation instruction [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/) \
Remember to append the path to the new BLAST bin directory to the existing PATH setting. 

## Usage

To run the script, go to the script folder and then use the following command:


```bash
python annotate_plasmids.py -i <input_file_or_directory> -o <output_directory> -t <file_type>
```

if you are running the script from other folder be sure to specify the path to annotate_plasmid.py. For example:

```bash
python /path_to_the_main_pipeline_folder/Scripts/annotate_plasmids.py -i <input_file_or_directory> -o <output_directory> -t <file_type>
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

The Databases can be downloaded from [here](https://rochester.box.com/v/PlasAnndatabases)
The first time the pipeline is run it will automatically download the databases.

Regarding any issues [Contact me](hislam2@ur.rochester.edu)
