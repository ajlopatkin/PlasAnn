
import sys
import os
import subprocess as sp
import time
import shutil
import argparse
import gdown

sys.path.append(os.path.join(os.path.dirname(__file__), 'Scripts'))
import essential_annotation
import draw_plasmid

def download_databases():
    # Get the absolute path to the Databases directory within the Scripts folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'Databases')
    
    # Check if the database directory already exists
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print("Database already exists. Skipping download.")
        return
    else:
        print("Downloading databases...")
        folder_id = '14jAiNrnsD7p0Kje--nB23fq_zTTE9noz'  # Folder ID extracted from your Google Drive link
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        gdown.download_folder(id=folder_id, output=output_dir, quiet=False)
        print("Download completed.")



def annotate_genbank_overwrite(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time = time.time()
    print(plasmid)
    CDS_list, DNA_seq, length_of_seq= essential_annotation.extract_genbank_info(plasmid,pathofdir)
    if DNA_seq is None:
        print(f"No sequence data available for {plasmid}. Skipping this file.")
        return 
    print(length_of_seq)
    essential_annotation.save_fasta_file(plasmid,DNA_seq)
    os.system(f'prodigal -i tmp_files/{plasmid}.fasta -o tmp_files/{plasmid}prodigal.gbk -f gbk > /dev/null 2>&1')
    positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
    complement_start, complement_end = essential_annotation.complementpositions(plasmid)
    DNA_CDS_list=essential_annotation.getlistofDNACDS(positions_of_coding_sequences,DNA_seq)
    essential_annotation.makequeryfastafordbsearch(DNA_CDS_list)
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    path_of_fasta='tmp_files'
    Initial_dataframe = essential_annotation.initial_blast_against_database(DNA_CDS_list, positions_of_coding_sequences, database)
    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, path_of_fasta)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, path_of_fasta)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, path_of_fasta)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, path_of_fasta, Initial_dataframe)
    final_dataframe = essential_annotation.process_final_dataframe(Initial_dataframe, replicon_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file(DNA_seq, final_dataframe, output_path_for_genbank, plasmid)
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path,plasmid)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")

# Similarly update annotate_genbank_retain and process_plasmid functions with the correct database paths

def annotate_genbank_retain(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time= time.time()
    print(plasmid)
    CDS_list, DNA_seq, length_of_seq= essential_annotation.extract_genbank_info(plasmid,pathofdir)
    if DNA_seq is None:
        print(f"No sequence data available for {plasmid}. Skipping this file.")
        return 
    DNA_CDS_list=essential_annotation.editing_cds_list(CDS_list)
    essential_annotation.save_fasta_file(plasmid,DNA_seq)
    path_of_fasta='tmp_files'
    positions_of_coding_sequences=essential_annotation.getpositionsofCDS_genbank(plasmid,pathofdir)
    complement_start,complement_end =essential_annotation.complementpositions_genbank(plasmid,pathofdir)
    Initial_dataframe = essential_annotation.initial_blast_against_database_genbank(DNA_CDS_list,positions_of_coding_sequences,default_db_path)
    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, path_of_fasta)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, path_of_fasta)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, path_of_fasta)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, path_of_fasta, Initial_dataframe)
    final_dataframe = essential_annotation.process_final_dataframe_genbank(Initial_dataframe, replicon_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file(DNA_seq, final_dataframe, output_path_for_genbank, plasmid)
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path,plasmid)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")


def process_plasmid(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time = time.time()
    print(plasmid)
    os.system(f'prodigal -i {os.path.join(pathofdir, plasmid)}.fasta -o tmp_files/{plasmid}prodigal.gbk -f gbk > /dev/null 2>&1')
    length_of_plasmid_sequence = essential_annotation.getsequencelength(plasmid, pathofdir)
    positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
    complement_start, complement_end = essential_annotation.complementpositions(plasmid)
    DNA_Sequence_of_fasta = essential_annotation.getthesequences(plasmid, pathofdir)
    DNA_cds_list = essential_annotation.getlistofDNACDS(positions_of_coding_sequences, DNA_Sequence_of_fasta)
    essential_annotation.makequeryfastafordbsearch(DNA_cds_list)
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    Initial_dataframe = essential_annotation.initial_blast_against_database(DNA_cds_list, positions_of_coding_sequences, database)
    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, pathofdir)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, pathofdir)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, pathofdir)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, pathofdir, Initial_dataframe)
    final_dataframe = essential_annotation.process_final_dataframe(Initial_dataframe, replicon_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file(DNA_Sequence_of_fasta, final_dataframe, output_path_for_genbank, plasmid)
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path,plasmid)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")

def main():
    parser = argparse.ArgumentParser(description='Annotate plasmid sequences from files.')
    parser.add_argument('-i', '--input', required=True, help='Input file or directory containing files.')
    parser.add_argument('-o', '--output', required=True, help='Output directory where the results will be stored.')
    parser.add_argument('-t', '--type', required=True, choices=['fasta', 'genbank'], help='Type of the input files either fasta or genbank.')

    args = parser.parse_args()

    # Download the databases automatically
    download_databases()

    # Get the absolute path to the Databases directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    databases_dir = os.path.join(script_dir, 'Databases')

    if args.type == 'genbank':
        choice = input("Choose an option:\n1. Retain existing CDS in GenBank files. This option wont use prodigal to detect CDS\n2. Overwrite existing CDS in GenBank files. This option will use prodigal to detect CDS\nEnter 1 or 2: ")
        if choice == '1':
            file_process_function = annotate_genbank_retain
        elif choice == '2':
            file_process_function = annotate_genbank_overwrite
        else:
            print("Invalid choice. Exiting...")
            sys.exit(1)
    else:
        file_process_function = process_plasmid

    if not os.path.exists("tmp_files"):
        os.makedirs("tmp_files")

    file_extension = '.fasta' if args.type == 'fasta' else '.gbk'

    if os.path.isdir(args.input):
        entries = os.listdir(args.input)
        file_list = [os.path.splitext(file)[0] for file in entries if os.path.isfile(os.path.join(args.input, file)) and file.endswith(file_extension)]
        for file_name in file_list:
            file_process_function(file_name, args.input, os.path.join(databases_dir, "DatabaseV5.1.3.csv"), 
                                  os.path.join(databases_dir, "oric.fna"), os.path.join(databases_dir, "orit.fna"), 
                                  os.path.join(databases_dir, "plasmidfinder.fasta"), os.path.join(databases_dir, "transposon.fasta"), args.output)
    elif os.path.isfile(args.input) and args.input.endswith(file_extension):
        file_name = os.path.splitext(os.path.basename(args.input))[0]
        file_process_function(file_name, os.path.dirname(args.input), os.path.join(databases_dir, "DatabaseV5.1.3.csv"), 
                              os.path.join(databases_dir, "oric.fna"), os.path.join(databases_dir, "orit.fna"), 
                              os.path.join(databases_dir, "plasmidfinder.fasta"), os.path.join(databases_dir, "transposon.fasta"), args.output)
    else:
        print("Invalid path or file type. Please provide a valid directory or file.")

    shutil.rmtree('tmp_files', ignore_errors=True)
    shutil.rmtree('makedb_folder', ignore_errors=True)

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
