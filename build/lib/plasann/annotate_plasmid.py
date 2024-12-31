
import os
import sys
import subprocess as sp
import time
import shutil
from . import essential_annotation
from . import draw_plasmid

def download_databases(output_dir):
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

def get_full_file_path(plasmid, pathofdir):
    # Define possible extensions
    extensions = ['.gbk', '.gb', '.genbank']
    for ext in extensions:
        file_path = os.path.join(pathofdir, f"{plasmid}{ext}")
        if os.path.exists(file_path):
            return file_path
    return None  # Return None if no file is found with any of the extensions

def annotate_genbank_overwrite(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time = time.time()
    print(plasmid)

    # Get the correct file path using the helper function
    file_path = get_full_file_path(plasmid, pathofdir)
    if file_path is None:
        print(f"No sequence data available for {plasmid}. Skipping this file.")
        return
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # Use the correct file path to extract GenBank information
    CDS_list, DNA_seq, length_of_seq = essential_annotation.extract_genbank_info(file_path)
    if DNA_seq =="":
        print(f"No sequence data available for {plasmid}.\n")
        print ("The pipeline will only annotate existing sequences in Genbank file.")
        print('No Transposable elements, ORIC, ORIT and Replicons will be annotated')
        print("To annotate these features use Fasta files as input or a genbank file with sequence content")
        updated_cds_list =essential_annotation.editing_cds_list(CDS_list)
        positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)
        complement_start, complement_end = essential_annotation.complementpositions_genbank(file_path,positions_of_coding_sequences)
        Initial_dataframe = essential_annotation.initial_blast_against_database_genbank(updated_cds_list, positions_of_coding_sequences, default_db_path)
        final_dataframe=essential_annotation.filter_close_sequences_cds(Initial_dataframe)
        #print(final_dataframe)
        final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
        output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
        essential_annotation.make_genbank_file_for_overwriting_cds_no_seq(length_of_seq, final_dataframe, output_path_for_genbank, plasmid,complement_start, complement_end)
        plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
        draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path, plasmid)
        essential_annotation.fix_genbank_date(output_path_for_genbank)
        end_time = time.time()
        duration = end_time - start_time
        print(f"The function took {duration} seconds to complete for {plasmid}.")
        return
    print(length_of_seq)
    essential_annotation.save_fasta_file(plasmid, DNA_seq)
    os.system(f'prodigal -i tmp_files/{plasmid}.fasta -o tmp_files/{plasmid}prodigal.gbk -f gbk -p meta > /dev/null 2>&1')
    positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
    complement_start, complement_end = essential_annotation.complementpositions(plasmid)
    DNA_CDS_list = essential_annotation.getlistofDNACDS(positions_of_coding_sequences, DNA_seq)
    essential_annotation.makequeryfastafordbsearch(DNA_CDS_list)
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    path_of_fasta = 'tmp_files'
    Initial_dataframe = essential_annotation.initial_blast_against_database(DNA_CDS_list, positions_of_coding_sequences, database)
    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, path_of_fasta)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, path_of_fasta)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, path_of_fasta)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, path_of_fasta, Initial_dataframe)
    final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(Initial_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe, replicon_dataframe)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file(DNA_seq, final_dataframe, output_path_for_genbank, plasmid)
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path, plasmid)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")



def annotate_genbank_retain(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time = time.time()
    print(plasmid)

    # Get the correct file path using the helper function
    file_path = get_full_file_path(plasmid, pathofdir)
    if file_path is None:
        print(f"No sequence data available for {plasmid}. Skipping this file.")
        return
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # Use the correct file path to extract GenBank information
    CDS_list, DNA_seq, length_of_seq = essential_annotation.extract_genbank_info(file_path)
    if DNA_seq =="":
        if DNA_seq =="":
            print(f"No sequence data available for {plasmid}.\n")
            print ("The pipeline will only annotate existing sequences in Genbank file.")
            print('No Transposable elements, ORIC, ORIT and Replicons will be annotated')
            print("To annotate these features use Fasta files as input or a genbank file with sequence content")
            updated_cds_list =essential_annotation.editing_cds_list(CDS_list)
            positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)
            complement_start, complement_end = essential_annotation.complementpositions_genbank(file_path,positions_of_coding_sequences)
            Initial_dataframe = essential_annotation.initial_blast_against_database_genbank(updated_cds_list, positions_of_coding_sequences, default_db_path)
            final_dataframe=essential_annotation.filter_close_sequences_cds(Initial_dataframe)
            #print(final_dataframe)
            final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
            output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
            essential_annotation.make_genbank_file_for_overwriting_cds_no_seq(length_of_seq, final_dataframe, output_path_for_genbank, plasmid,complement_start, complement_end)
            plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
            draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path, plasmid)
            essential_annotation.fix_genbank_date(output_path_for_genbank)
            end_time = time.time()
            duration = end_time - start_time
            print(f"The function took {duration} seconds to complete for {plasmid}.")
        return
    DNA_CDS_list = essential_annotation.editing_cds_list(CDS_list)
    essential_annotation.save_fasta_file(plasmid, DNA_seq)
    path_of_fasta = 'tmp_files'
    positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)
    complement_start, complement_end = essential_annotation.complementpositions_genbank(file_path,positions_of_coding_sequences)
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    Initial_dataframe = essential_annotation.initial_blast_against_database_genbank(DNA_CDS_list, positions_of_coding_sequences, default_db_path)
    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, path_of_fasta)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, path_of_fasta)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, path_of_fasta)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, path_of_fasta, Initial_dataframe)

    final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(Initial_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe, replicon_dataframe)
    final_dataframe.to_csv('/Users/habibulislam/Desktop/Database_stats/test/F_final.csv', index=False)
    output_path = os.path.join(output_directory, plasmid)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file_for_retaining_cds(DNA_seq, final_dataframe, output_path_for_genbank, plasmid,complement_start, complement_end)
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path, plasmid)
    essential_annotation.fix_genbank_date(output_path_for_genbank)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")


def process_plasmid(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    start_time = time.time()
    print("Starting the annotation process of " + plasmid)

    # List of possible FASTA file extensions
    extensions = ['.fasta', '.fa', '.fsa', '.fna']
    fasta_file_path = None

    # Try to find the file with one of the specified extensions
    for ext in extensions:
        potential_path = os.path.join(pathofdir, plasmid + ext)
        if os.path.exists(potential_path):
            fasta_file_path = potential_path
            print(fasta_file_path)
            break

    # Check if a valid FASTA file was found
    '''if not fasta_file_path:
        print(f"No valid FASTA file found for {plasmid} with any known extensions.")
        return  # Exit the function early if no valid FASTA file is found'''

    os.system(f'prodigal -i "{fasta_file_path}" -o "tmp_files/{plasmid}prodigal.gbk" -f gbk -p meta > /dev/null 2>&1')

    length_of_plasmid_sequence = essential_annotation.getsequencelength(plasmid, pathofdir)
    print(f"Length of plasmid {plasmid}: {length_of_plasmid_sequence}")
    gc_content = essential_annotation.calculate_gc_content(fasta_file_path)
    print(f"GC Content of Plasmid {plasmid}: {gc_content:.2f}%")

    positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
    complement_start, complement_end = essential_annotation.complementpositions(plasmid)
    DNA_Sequence_of_fasta = essential_annotation.getthesequences(plasmid, pathofdir)
    DNA_cds_list = essential_annotation.getlistofDNACDS(positions_of_coding_sequences, DNA_Sequence_of_fasta)

    if not DNA_cds_list:
        print("Prodigal could not detect any coding sequences.")
        return  # Exit the function early if no CDS are found

    essential_annotation.makequeryfastafordbsearch(DNA_cds_list)
    database = essential_annotation.makedatabasefromcsvfile(default_db_path)
    Initial_dataframe = essential_annotation.initial_blast_against_database(DNA_cds_list, positions_of_coding_sequences, database)

    oric_dataframe = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, pathofdir)
    orit_dataframe = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, pathofdir)
    replicon_dataframe = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, pathofdir)
    transposon_dataframe = essential_annotation.blast_against_transposon_database(transposon_data, plasmid, pathofdir, Initial_dataframe)

    final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(Initial_dataframe, oric_dataframe, orit_dataframe, transposon_dataframe, replicon_dataframe)
    output_path = os.path.join(output_directory, plasmid)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    final_dataframe.to_csv(os.path.join(output_path, f"Annotation_table_for_{plasmid}.csv"))
    output_path_for_genbank = os.path.join(output_path, f"Annotation_gbk_file_for_{plasmid}.gbk")
    essential_annotation.make_genbank_file(DNA_Sequence_of_fasta, final_dataframe, output_path_for_genbank, plasmid)
    
    plasmid_map_path = os.path.join(output_path, f"Annotated_Map_for_{plasmid}.png")
    draw_plasmid.draw_plasmid_map_from_genbank_file(output_path_for_genbank, plasmid_map_path, plasmid)
    essential_annotation.fix_genbank_date(output_path_for_genbank)
    end_time = time.time()
    duration = end_time - start_time
    print(f"The function took {duration} seconds to complete for {plasmid}.")


