# src/plasann/cli.py
import argparse
import os
import shutil
import sys
from pathlib import Path
from . import annotate_plasmid  # This imports the module
import gdown

def download_databases(output_dir):
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print("Database already exists. Skipping download.")
        return
    else:
        print("Downloading databases...")
        folder_id = '14jAiNrnsD7p0Kje--nB23fq_zTTE9noz'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        gdown.download_folder(id=folder_id, output=output_dir, quiet=False)
        print("Download completed.")

def main():
    parser = argparse.ArgumentParser(description='Annotate plasmid sequences from files.')
    parser.add_argument('-i', '--input', required=True, 
                       help='Input file or directory containing files.')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output directory where the results will be stored.')
    parser.add_argument('-t', '--type', required=True, 
                       choices=['fasta', 'genbank'], 
                       help='Type of the input files either fasta or genbank.')

    args = parser.parse_args()

    # Download the databases automatically
    databases_dir = "Databases"
    download_databases(databases_dir)

    if args.type == 'genbank':
        choice = input("Choose an option:\n1. Retain existing CDS in GenBank files\n2. Overwrite existing CDS in GenBank files\nEnter 1 or 2: ")
        if choice == '1':
            file_process_function = annotate_plasmid.annotate_genbank_retain  # Changed this line
        elif choice == '2':
            file_process_function = annotate_plasmid.annotate_genbank_overwrite  # Changed this line
        else:
            print("Invalid choice. Exiting...")
            sys.exit(1)
    else:
        file_process_function = annotate_plasmid.process_plasmid  # Changed this line

    # Create tmp_files directory if it doesn't exist
    if not os.path.exists("tmp_files"):
        os.makedirs("tmp_files")

    # Rest of your code remains the same
    # Just make sure all function calls use annotate_plasmid.function_name
    if os.path.isdir(args.input):
        entries = os.listdir(args.input)
        file_list = [os.path.join(args.input, file) for file in entries 
                    if os.path.isfile(os.path.join(args.input, file))]
        for file_path in file_list:
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            file_process_function(file_name, os.path.dirname(file_path),
                                os.path.join(databases_dir, "Database.csv"),
                                os.path.join(databases_dir, "oric.fna"),
                                os.path.join(databases_dir, "orit.fna"),
                                os.path.join(databases_dir, "plasmidfinder.fasta"),
                                os.path.join(databases_dir, "transposon.fasta"),
                                args.output)
    elif os.path.isfile(args.input):
        file_name = os.path.splitext(os.path.basename(args.input))[0]
        file_process_function(file_name, os.path.dirname(args.input),
                            os.path.join(databases_dir, "Database.csv"),
                            os.path.join(databases_dir, "oric.fna"),
                            os.path.join(databases_dir, "orit.fna"),
                            os.path.join(databases_dir, "plasmidfinder.fasta"),
                            os.path.join(databases_dir, "transposon.fasta"),
                            args.output)
    else:
        print("Invalid path or file type. Please provide a valid directory or file.")

    # Clean up
    shutil.rmtree('tmp_files', ignore_errors=True)
    shutil.rmtree('makedb_folder', ignore_errors=True)

if __name__ == "__main__":
    main()