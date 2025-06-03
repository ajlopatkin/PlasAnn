# PlasAnn

A comprehensive plasmid annotation tool that enables the identification and visualization of key plasmid features.

## Installation
pip install plasann

## Usage
PlasAnn -i input.fasta -o output_directory -t fasta

or for GenBank files:
PlasAnn -i input.gbk -o output_directory -t genbank

## Features

- Automated identification of coding sequences
- Detection of origins of replication (oriC)
- Detection of origins of transfer (oriT)
- Identification of mobile genetic elements
- Recognition of replicon types
- Generation of annotated GenBank files
- Creation of circular plasmid maps

## Dependencies

- BLAST+ (must be installed separately)
- Prodigal (must be installed separately)

## License

MIT