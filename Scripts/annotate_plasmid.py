import sys
import os
import subprocess as sp
import time
import shutil
import argparse
import gdown

from . import essential_annotation
from . import draw_plasmid
sys.path.append(os.path.join(os.path.dirname(__file__), 'Scripts'))

import warnings
from Bio import BiopythonParserWarning

# Terminal colors for better visibility
class TermColors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# Suppress Biopython parser warnings
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

def print_header(message):
    """Print a header message with decoration"""
    print("\n" + "=" * 80)
    print(f"{TermColors.BOLD}{TermColors.HEADER}{message}{TermColors.ENDC}")
    print("=" * 80)

def print_subheader(message):
    """Print a subheader message with decoration"""
    print(f"\n{TermColors.BOLD}{TermColors.BLUE}>> {message}{TermColors.ENDC}")
    print("-" * 60)

def print_success(message):
    """Print a success message"""
    print(f"{TermColors.GREEN}✓ {message}{TermColors.ENDC}")

def print_warning(message):
    """Print a warning message"""
    print(f"{TermColors.YELLOW}⚠ {message}{TermColors.ENDC}")

def print_error(message):
    """Print an error message"""
    print(f"{TermColors.RED}✗ {message}{TermColors.ENDC}")

def print_info(message):
    """Print an info message"""
    print(f"{TermColors.CYAN}ℹ {message}{TermColors.ENDC}")

def print_step(step_number, message):
    """Print a step in a process"""
    print(f"{TermColors.BOLD}[{step_number}] {message}{TermColors.ENDC}")

def print_stats(label, value, unit=""):
    """Print a statistic with label and value"""
    print(f"  {TermColors.BOLD}{label}:{TermColors.ENDC} {value}{unit}")



def find_fasta_file(plasmid, directory):
    """
    Find a FASTA file for a given plasmid name in the specified directory.
    
    Args:
        plasmid (str): The name of the plasmid (without extension)
        directory (str): The directory to search in
        
    Returns:
        str or None: Path to the FASTA file if found, None otherwise
    """
    extensions = ['.fasta', '.fa', '.fsa', '.fna']
    
    for ext in extensions:
        potential_path = os.path.join(directory, plasmid + ext)
        if os.path.exists(potential_path):
            return potential_path
    
    return None


def run_prodigal(input_file, output_file):
    """
    Run Prodigal to predict coding sequences in the input file.
    
    Args:
        input_file (str): Path to the input FASTA file
        output_file (str): Path to the output GenBank file
        
    Returns:
        int: Return code from Prodigal command
    """
    command = f'prodigal -i "{input_file}" -o "{output_file}" -f gbk -p meta > /dev/null 2>&1'
    return os.system(command)


def download_databases(output_dir):
    # Check if the database directory already exists
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print_success(f"Database already exists in {output_dir}. Skipping download.")
        return
    else:
        print_header("Downloading Annotation Databases")
        print_info("This may take a few minutes depending on your internet connection...")
        
        folder_id = '14jAiNrnsD7p0Kje--nB23fq_zTTE9noz'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        gdown.download_folder(id=folder_id, output=output_dir, quiet=False)
        print_success("Database download completed successfully!")

def process_plasmid(plasmid, pathofdir, default_db_path, oric_data, orit_data, plasmidfinder_data, transposon_data, output_directory):
    """
    Process a plasmid from a FASTA file, extract coding sequences, and annotate the plasmid.
    """
    start_time = time.time()
    print_header(f"Annotating Plasmid: {plasmid}")
    
    # Find the FASTA file
    print_step(1, "Locating FASTA file")
    fasta_file_path = find_fasta_file(plasmid, pathofdir)
    if not fasta_file_path:
        print_error(f"No valid FASTA file found for {plasmid} with any known extensions.")
        return
    print_success(f"Found: {fasta_file_path}")
    
    # Validate file isn't empty or corrupted
    try:
        file_size = os.path.getsize(fasta_file_path)
        if file_size == 0:
            print_error(f"Error: {fasta_file_path} is empty. Skipping this file.")
            return
            
        # Quickly check if file is valid FASTA format
        with open(fasta_file_path, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                print_error(f"Error: {fasta_file_path} doesn't appear to be a valid FASTA file. Skipping.")
                return
    except Exception as e:
        print_error(f"Error reading {fasta_file_path}: {e}")
        return
        
    # Create output directory for this plasmid
    plasmid_output_dir = os.path.join(output_directory, plasmid)
    if not os.path.exists(plasmid_output_dir):
        os.makedirs(plasmid_output_dir)
    
    # Run prodigal to identify coding sequences
    print_step(2, "Running Prodigal for gene prediction")
    prodigal_output = f"tmp_files/{plasmid}prodigal.gbk"
    run_result = run_prodigal(fasta_file_path, prodigal_output)
    
    if run_result != 0:
        print_error(f"Prodigal failed to run for {plasmid}. Return code: {run_result}")
        return
    print_success("Prodigal completed successfully")
    
    # Get basic plasmid information
    print_step(3, "Analyzing plasmid sequence properties")
    try:
        dna_sequence = essential_annotation.getthesequences(plasmid, pathofdir)
        if dna_sequence == "File not found":
            print_error(f"Error retrieving sequence for {plasmid}. Aborting.")
            return
            
        sequence_length = len(dna_sequence)
        gc_content = essential_annotation.calculate_gc_content(fasta_file_path)
        
        print_stats("Length", sequence_length, " bp")
        print_stats("GC Content", f"{gc_content:.2f}", "%")
    except Exception as e:
        print_error(f"Error processing sequence data for {plasmid}: {e}")
        return
    
    # Extract coding sequences
    print_step(4, "Extracting coding sequences")
    try:
        positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
        complement_start, complement_end = essential_annotation.complementpositions(plasmid)
        dna_cds_list = essential_annotation.getlistofDNACDS(positions_of_coding_sequences, dna_sequence)
        print_info(f"Found {len(dna_cds_list)} potential coding sequences")
    except Exception as e:
        print_error(f"Error extracting coding sequences for {plasmid}: {e}")
        return
    
    if not dna_cds_list:
        print_warning(f"Prodigal could not detect any coding sequences in {plasmid}.")
        return
    
    # Prepare database and search
    print_step(5, "Performing sequence annotation against databases")
    try:
        print_info("Preparing query sequences...")
        essential_annotation.makequeryfastafordbsearch(dna_cds_list)
        database = essential_annotation.makedatabasefromcsvfile(default_db_path)
        
        print_info("Searching against main annotation database...")
        initial_dataframe = essential_annotation.initial_blast_against_database(
            dna_cds_list, positions_of_coding_sequences, database)
        #initial_dataframe.to_csv('initial_dataframe.csv')
        
        # Annotate with special features
        print_info("Searching for origin of replication (OriC)...")
        oric_df = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, pathofdir)
        #oric_df.to_csv('oric.csv')
        
        print_info("Searching for origin of transfer (OriT)...")
        orit_df = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, pathofdir)
        #orit_df.to_csv('orit.csv')
        
        print_info("Searching for replicon types...")
        replicon_df = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, pathofdir)
        #replicon_df.to_csv('replicon.csv')
        
        print_info("Searching for mobile genetic elements...")
        transposon_df = essential_annotation.blast_against_transposon_database(
            transposon_data, plasmid, pathofdir, initial_dataframe)
        #transposon_df.to_csv('transposon.csv')
    except Exception as e:
        print_error(f"Error during database search for {plasmid}: {e}")
        return
    
    # Merge all annotations
    print_step(6, "Finalizing annotations and generating output files")
    try:
        print_info("Merging all annotation data...")
        final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(
            initial_dataframe, oric_df, orit_df, transposon_df, replicon_df)
        
        # Save annotation results
        annotation_table_path = os.path.join(plasmid_output_dir, f"Annotation_table_for_{plasmid}.csv")
        final_dataframe.to_csv(annotation_table_path)
        print_success(f"Saved annotation table to: {annotation_table_path}")
        
        # Create GenBank file
        print_info("Creating GenBank file...")
        genbank_output_path = os.path.join(plasmid_output_dir, f"Annotation_gbk_file_for_{plasmid}.gbk")
        essential_annotation.make_genbank_file(dna_sequence, final_dataframe, genbank_output_path, plasmid)
        essential_annotation.adjust_start_positions_in_place(genbank_output_path)
        print_success(f"Saved GenBank file to: {genbank_output_path}")

        
        # Create plasmid map
        print_info("Generating plasmid map...")
        plasmid_map_path = os.path.join(plasmid_output_dir, f"Annotated_Map_for_{plasmid}.png")
        draw_plasmid.draw_plasmid_map_from_genbank_file(genbank_output_path, plasmid_map_path, plasmid)
        print_success(f"Saved plasmid map to: {plasmid_map_path}")
        
        # Fix the date in the GenBank file
        essential_annotation.fix_genbank_date(genbank_output_path)
    except Exception as e:
        print_error(f"Error finalizing annotation for {plasmid}: {e}")
        return
    
    # Report execution time
    end_time = time.time()
    duration = end_time - start_time
    print_header(f"Annotation Complete: {plasmid}")
    print_stats("Processing time", f"{duration:.2f}", " seconds")
    print_success(f"All files saved to: {plasmid_output_dir}")

def get_full_file_path(plasmid, pathofdir):
    """
    Find a GenBank file for a given plasmid name in the specified directory.
    
    Args:
        plasmid (str): The name of the plasmid (without extension)
        directory (str): The directory to search in
        
    Returns:
        str or None: Path to the GenBank file if found, None otherwise
    """
    extensions = ['.gbk', '.gb', '.genbank']
    
    for ext in extensions:
        potential_path = os.path.join(pathofdir, plasmid + ext)
        if os.path.exists(potential_path):
            return potential_path
    
    return None

def annotate_genbank_overwrite(plasmid, pathofdir, default_db_path, oric_data, orit_data, 
                               plasmidfinder_data, transposon_data, output_directory):
    """
    Process a GenBank file by extracting its sequence and overwriting CDS annotations.
    If no sequence is available, tries to fetch from NCBI or annotate existing CDS.
    
    Args:
        plasmid (str): The name of the plasmid (without extension)
        pathofdir (str): Path to the directory containing the plasmid file
        default_db_path (str): Path to the default database
        oric_data (str): Path to the OriC database
        orit_data (str): Path to the OriT database
        plasmidfinder_data (str): Path to the PlasmidFinder database
        transposon_data (str): Path to the transposon database
        output_directory (str): Path to the output directory
    """
    start_time = time.time()
    print_header(f"Annotating GenBank: {plasmid} (Overwrite Mode)")
    
    # Find the GenBank file
    print_step(1, "Locating GenBank file")
    file_path = get_full_file_path(plasmid, pathofdir)
    if file_path is None:
        print_error(f"No GenBank file found for {plasmid}. Skipping this file.")
        return
    print_success(f"Found: {file_path}")
    
    # Create output directory for this plasmid
    plasmid_output_dir = os.path.join(output_directory, plasmid)
    if not os.path.exists(plasmid_output_dir):
        os.makedirs(plasmid_output_dir)
    
    # Extract GenBank information
    print_step(2, "Extracting sequence information from GenBank")
    try:
        CDS_list, DNA_seq, length_of_seq = essential_annotation.extract_genbank_info(file_path)
        print_stats("GenBank file length", length_of_seq, " bp")
    except Exception as e:
        print_error(f"Error extracting information from GenBank file: {e}")
        return
    
    # Check if sequence is available
    if not DNA_seq:
        print_warning("No sequence data available in the GenBank file.")
        print_info("Attempting to process existing GenBank annotations...")
         
        try:
            # Process existing annotations when no sequence is available
            print_step(3, "Processing existing annotations")
            CDS_info= essential_annotation.extract_genbank_edited(file_path)
            list_of_cds = [cds[2] for cds in CDS_info if cds[2]]  # Get sequences
            list_of_positions = [(cds[0][0], cds[0][1]) for cds in CDS_info if cds[2]]
            print_info(f"Found {len(list_of_cds)} CDS features in GenBank file")
            #print(CDS_info)
            

            '''updated_cds_list = essential_annotation.editing_cds_list(CDS_list)
            positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)'''
            
            print_step(4, "Identifying complement regions")
            complement_start, complement_end = essential_annotation.extract_reverse_positions(CDS_info)
            
            
            print_step(5, "Annotating CDS regions")
            database = essential_annotation.makedatabasefromcsvfile(default_db_path)
            initial_dataframe = essential_annotation.initial_blast_against_database_genbank(
                list_of_cds, list_of_positions, default_db_path)
            final_dataframe = essential_annotation.filter_close_sequences_cds(initial_dataframe)
            
            # Save annotation results
            print_step(6, "Saving annotation results")
            annotation_table_path = os.path.join(plasmid_output_dir, f"Annotation_table_for_{plasmid}.csv")
            final_dataframe.to_csv(annotation_table_path)
            print_success(f"Saved annotation table to: {annotation_table_path}")
            
            # Create GenBank file
            print_info("Creating GenBank file...")
            genbank_output_path = os.path.join(plasmid_output_dir, f"Annotation_gbk_file_for_{plasmid}.gbk")
            essential_annotation.make_genbank_file_for_overwriting_cds_no_seq(
                length_of_seq, final_dataframe, genbank_output_path, plasmid, complement_start, complement_end)
            print_success(f"Saved GenBank file to: {genbank_output_path}")
            
            # Create plasmid map
            print_info("Generating plasmid map...")
            plasmid_map_path = os.path.join(plasmid_output_dir, f"Annotated_Map_for_{plasmid}.png")
            draw_plasmid.draw_plasmid_map_from_genbank_file(genbank_output_path, plasmid_map_path, plasmid)
            print_success(f"Saved plasmid map to: {plasmid_map_path}")
            
            # Fix GenBank date
            essential_annotation.fix_genbank_date(genbank_output_path)
            
        except Exception as e:
            print_error(f"Error processing existing annotations: {e}")
            return
    else:
        print_success("Sequence found in GenBank file. Proceeding with full annotation.")
        
        try:
            # Process with full sequence
            print_step(3, "Saving sequence as FASTA")
            essential_annotation.save_fasta_file(plasmid, DNA_seq)
            
            # Calculate GC content using the temporary FASTA file
            temp_fasta_path = f"tmp_files/{plasmid}.fasta"
            gc_percentage = essential_annotation.calculate_gc_content(temp_fasta_path)
            print_stats("GC Content", f"{gc_percentage:.2f}", "%")
            
            print_step(4, "Running Prodigal for gene prediction")
            run_result = run_prodigal(f"tmp_files/{plasmid}.fasta", f"tmp_files/{plasmid}prodigal.gbk")
            if run_result != 0:
                print_error(f"Prodigal failed to run for {plasmid}. Return code: {run_result}")
                return
            print_success("Prodigal completed successfully")
            
            print_step(5, "Extracting coding sequences")
            positions_of_coding_sequences = essential_annotation.getpositionsofCDS(plasmid)
            complement_start, complement_end = essential_annotation.complementpositions(plasmid)
            DNA_CDS_list = essential_annotation.getlistofDNACDS(positions_of_coding_sequences, DNA_seq)
            print_info(f"Found {len(DNA_CDS_list)} potential coding sequences")
            
            print_step(6, "Performing sequence annotation")
            essential_annotation.makequeryfastafordbsearch(DNA_CDS_list)
            database = essential_annotation.makedatabasefromcsvfile(default_db_path)
            
            print_info("Searching against main annotation database...")
            initial_dataframe = essential_annotation.initial_blast_against_database(
                DNA_CDS_list, positions_of_coding_sequences, database)
            
            # Annotate with special features
            print_info("Searching for origin of replication (OriC)...")
            oric_df = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, 'tmp_files')
            
            print_info("Searching for origin of transfer (OriT)...")
            orit_df = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, 'tmp_files')
            
            print_info("Searching for replicon types...")
            replicon_df = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, 'tmp_files')
            
            print_info("Searching for mobile genetic elements...")
            transposon_df = essential_annotation.blast_against_transposon_database(
                transposon_data, plasmid, 'tmp_files', initial_dataframe)
            
            print_step(7, "Merging annotation data")
            final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(
                initial_dataframe, oric_df, orit_df, transposon_df, replicon_df)
            
            # Save annotation results
            print_step(8, "Saving results")
            annotation_table_path = os.path.join(plasmid_output_dir, f"Annotation_table_for_{plasmid}.csv")
            final_dataframe.to_csv(annotation_table_path)
            print_success(f"Saved annotation table to: {annotation_table_path}")
            
            # Create GenBank file
            print_info("Creating GenBank file...")
            genbank_output_path = os.path.join(plasmid_output_dir, f"Annotation_gbk_file_for_{plasmid}.gbk")
            essential_annotation.make_genbank_file(DNA_seq, final_dataframe, genbank_output_path, plasmid)
            essential_annotation.adjust_start_positions_in_place(genbank_output_path)
            print_success(f"Saved GenBank file to: {genbank_output_path}")
            
            # Create plasmid map
            print_info("Generating plasmid map...")
            plasmid_map_path = os.path.join(plasmid_output_dir, f"Annotated_Map_for_{plasmid}.png")
            draw_plasmid.draw_plasmid_map_from_genbank_file(genbank_output_path, plasmid_map_path, plasmid)
            print_success(f"Saved plasmid map to: {plasmid_map_path}")
            
            # Fix GenBank date
            essential_annotation.fix_genbank_date(genbank_output_path)
            
        except Exception as e:
            print_error(f"Error during annotation process: {e}")
            return
    
    # Report execution time
    end_time = time.time()
    duration = end_time - start_time
    print_header(f"Annotation Complete: {plasmid}")
    print_stats("Processing time", f"{duration:.2f}", " seconds")
    print_success(f"All files saved to: {plasmid_output_dir}")

def annotate_genbank_retain(plasmid, pathofdir, default_db_path, oric_data, orit_data, 
                            plasmidfinder_data, transposon_data, output_directory):
    """
    Process a GenBank file by extracting its sequence and retaining original CDS annotations.
    If no sequence is available, tries to fetch from NCBI or annotate existing CDS.
    
    Args:
        plasmid (str): The name of the plasmid (without extension)
        pathofdir (str): Path to the directory containing the plasmid file
        default_db_path (str): Path to the default database
        oric_data (str): Path to the OriC database
        orit_data (str): Path to the OriT database
        plasmidfinder_data (str): Path to the PlasmidFinder database
        transposon_data (str): Path to the transposon database
        output_directory (str): Path to the output directory
    """
    start_time = time.time()
    print_header(f"Annotating GenBank: {plasmid} (Retain Mode)")
    
    # Find the GenBank file
    print_step(1, "Locating GenBank file")
    file_path = get_full_file_path(plasmid, pathofdir)
    if file_path is None:
        print_error(f"No GenBank file found for {plasmid}. Skipping this file.")
        return
    print_success(f"Found: {file_path}")
    
    # Create output directory for this plasmid
    plasmid_output_dir = os.path.join(output_directory, plasmid)
    if not os.path.exists(plasmid_output_dir):
        os.makedirs(plasmid_output_dir)
    
    # Extract GenBank information
    print_step(2, "Extracting sequence information from GenBank")
    try:
        CDS_list, DNA_seq, length_of_seq = essential_annotation.extract_genbank_info(file_path)
        CDS_info = essential_annotation.extract_genbank_edited(file_path)
        print_stats("GenBank file length", length_of_seq, " bp")
    except Exception as e:
        print_error(f"Error extracting information from GenBank file: {e}")
        return
    
    # Check if sequence is available
    if not DNA_seq:
        print_warning("No sequence data available in the GenBank file.")
        print_info("Attempting to process existing GenBank annotations...")
        
        try:
            # Process existing annotations when no sequence is available
            print_step(3, "Processing existing annotations")
            CDS_info= essential_annotation.extract_genbank_edited(file_path)
            list_of_cds = [cds[2] for cds in CDS_info if cds[2]]  # Get sequences
            list_of_positions = [(cds[0][0], cds[0][1]) for cds in CDS_info if cds[2]]
            print_info(f"Found {len(list_of_cds)} CDS features in GenBank file")
            #print(CDS_info)
            

            '''updated_cds_list = essential_annotation.editing_cds_list(CDS_list)
            positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)'''
            
            print_step(4, "Identifying complement regions")
            complement_start, complement_end = essential_annotation.extract_reverse_positions(CDS_info)
            
            
            print_step(5, "Annotating CDS regions")
            database = essential_annotation.makedatabasefromcsvfile(default_db_path)
            initial_dataframe = essential_annotation.initial_blast_against_database_genbank(
                list_of_cds, list_of_positions, default_db_path)
            final_dataframe = essential_annotation.filter_close_sequences_cds(initial_dataframe)
            
            # Save annotation results
            print_step(6, "Saving annotation results")
            annotation_table_path = os.path.join(plasmid_output_dir, f"Annotation_table_for_{plasmid}.csv")
            final_dataframe.to_csv(annotation_table_path)
            print_success(f"Saved annotation table to: {annotation_table_path}")
            
            # Create GenBank file
            print_info("Creating GenBank file...")
            genbank_output_path = os.path.join(plasmid_output_dir, f"Annotation_gbk_file_for_{plasmid}.gbk")
            essential_annotation.make_genbank_file_for_overwriting_cds_no_seq(
                length_of_seq, final_dataframe, genbank_output_path, plasmid, complement_start, complement_end)
            print_success(f"Saved GenBank file to: {genbank_output_path}")
            
            # Create plasmid map
            print_info("Generating plasmid map...")
            plasmid_map_path = os.path.join(plasmid_output_dir, f"Annotated_Map_for_{plasmid}.png")
            draw_plasmid.draw_plasmid_map_from_genbank_file(genbank_output_path, plasmid_map_path, plasmid)
            print_success(f"Saved plasmid map to: {plasmid_map_path}")
            
            # Fix GenBank date
            essential_annotation.fix_genbank_date(genbank_output_path)
            
        except Exception as e:
            print_error(f"Error processing existing annotations: {e}")
            return
    else:
        print_success("Sequence found in GenBank file. Proceeding with full annotation.")
        
        try:
            # Process with full sequence and retain CDS
            print_step(3, "Extracting CDS information")
            list_of_cds = [cds[2] for cds in CDS_info if cds[2]]  # Get sequences
            list_of_positions = [(cds[0][0], cds[0][1]) for cds in CDS_info if cds[2]]  # Get positions
            
            print_info(f"Found {len(list_of_cds)} CDS features in GenBank file")
            
            print_step(4, "Saving sequence as FASTA")
            essential_annotation.save_fasta_file(plasmid, DNA_seq)
            
            print_step(5, "Identifying complement regions")
            positions_of_coding_sequences = essential_annotation.getpositionsofCDS_genbank(file_path)
            complement_start, complement_end = essential_annotation.complementpositions_genbank(
                file_path, list_of_positions)
            
            print_step(6, "Performing sequence annotation")
            database = essential_annotation.makedatabasefromcsvfile(default_db_path)
            initial_dataframe = essential_annotation.initial_blast_against_database_genbank(
                list_of_cds, list_of_positions, default_db_path)
            
            # Annotate with special features
            print_info("Searching for origin of replication (OriC)...")
            oric_df = essential_annotation.blast_against_oric_dataframe(oric_data, plasmid, 'tmp_files')
            
            print_info("Searching for origin of transfer (OriT)...")
            orit_df = essential_annotation.blast_against_orit_dataframe(orit_data, plasmid, 'tmp_files')
            
            print_info("Searching for replicon types...")
            replicon_df = essential_annotation.blast_against_replicon_database(plasmidfinder_data, plasmid, 'tmp_files')
            
            print_info("Searching for mobile genetic elements...")
            transposon_df = essential_annotation.blast_against_transposon_database(
                transposon_data, plasmid, 'tmp_files', initial_dataframe)
            
            print_step(7, "Merging annotation data")
            final_dataframe = essential_annotation.merge_all_ther_database_and_fix_accordingly(
                initial_dataframe, oric_df, orit_df, transposon_df, replicon_df)
            
            # Save annotation results
            print_step(8, "Saving results")
            annotation_table_path = os.path.join(plasmid_output_dir, f"Annotation_table_for_{plasmid}.csv")
            final_dataframe.to_csv(annotation_table_path)
            print_success(f"Saved annotation table to: {annotation_table_path}")
            
            # Create GenBank file
            print_info("Creating GenBank file...")
            genbank_output_path = os.path.join(plasmid_output_dir, f"Annotation_gbk_file_for_{plasmid}.gbk")
            essential_annotation.make_genbank_file_for_retaining_cds(
                DNA_seq, final_dataframe, genbank_output_path, plasmid, complement_start, complement_end)
            essential_annotation.update_genbank_file_with_reverse(genbank_output_path, CDS_info)
            essential_annotation.adjust_start_positions_in_place(genbank_output_path)
            print_success(f"Saved GenBank file to: {genbank_output_path}")
            
            # Create plasmid map
            print_info("Generating plasmid map...")
            plasmid_map_path = os.path.join(plasmid_output_dir, f"Annotated_Map_for_{plasmid}.png")
            draw_plasmid.draw_plasmid_map_from_genbank_file(genbank_output_path, plasmid_map_path, plasmid)
            print_success(f"Saved plasmid map to: {plasmid_map_path}")
            
            # Fix GenBank date
            essential_annotation.fix_genbank_date(genbank_output_path)
            
        except Exception as e:
            print_error(f"Error during annotation process: {e}")
            return
    
    # Report execution time
    end_time = time.time()
    duration = end_time - start_time
    print_header(f"Annotation Complete: {plasmid}")
    print_stats("Processing time", f"{duration:.2f}", " seconds")
    print_success(f"All files saved to: {plasmid_output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Annotate plasmid sequences from files.')
    parser.add_argument('-i', '--input', required=True, help='Input file or directory containing files.')
    parser.add_argument('-o', '--output', required=True, help='Output directory where the results will be stored.')
    parser.add_argument('-t', '--type', required=True, choices=['fasta', 'genbank'], help='Type of the input files either fasta or genbank.')

    args = parser.parse_args()
    
    print_header("Plasmid Annotation Pipeline")
    print_info(f"Input: {args.input}")
    print_info(f"Output: {args.output}")
    print_info(f"File type: {args.type}")

    # Download the databases automatically
    databases_dir = "Databases"
    download_databases(databases_dir)

    if args.type == 'genbank':
        print_subheader("GenBank Processing Options")
        choice = input("Choose an option:\n1. Retain existing CDS in GenBank files. This option won't use prodigal to detect CDS\n2. Overwrite existing CDS in GenBank files. This option will use prodigal to detect CDS\nEnter 1 or 2: ")
        if choice == '1':
            file_process_function = annotate_genbank_retain
            print_info("Selected: Retain existing CDS annotations")
        elif choice == '2':
            file_process_function = annotate_genbank_overwrite
            print_info("Selected: Overwrite existing CDS annotations using Prodigal")
        else:
            print_error("Invalid choice. Exiting...")
            sys.exit(1)
    else:
        file_process_function = process_plasmid
        print_info("Processing FASTA files using Prodigal for gene prediction")

    # Create temp directory if it doesn't exist
    if not os.path.exists("tmp_files"):
        os.makedirs("tmp_files")
        print_info("Created temporary directory: tmp_files")

    # Define the valid extensions for GenBank files
    genbank_extensions = ['.gbk', '.gb', '.genbank']

    # Define the valid extensions for Fasta files
    fasta_extensions = ['.fasta', '.fa', '.fsa', '.fna']

    # Determine the correct extension set based on file type
    valid_extensions = genbank_extensions if args.type == 'genbank' else fasta_extensions

    # Process input files
    if os.path.isdir(args.input):
        entries = os.listdir(args.input)
        # Explicitly filter out .DS_Store files and only include valid file extensions
        file_list = []
        for file in entries:
            # Skip .DS_Store files and any other hidden files
            if file.startswith('.'):
                continue
                
            file_path = os.path.join(args.input, file)
            if os.path.isfile(file_path) and any(file.endswith(ext) for ext in valid_extensions):
                file_list.append(file_path)
        
        print_subheader(f"Found {len(file_list)} valid files to process")
        
        for i, file_path in enumerate(file_list):
            file_name = os.path.splitext(os.path.basename(file_path))[0]
            print_subheader(f"Processing file {i+1} of {len(file_list)}: {file_name}")
            try:
                file_process_function(file_name, os.path.dirname(file_path), os.path.join(databases_dir, "Database.csv"),
                                    os.path.join(databases_dir, "oric.fna"), os.path.join(databases_dir, "orit.fna"),
                                    os.path.join(databases_dir, "plasmidfinder.fasta"), os.path.join(databases_dir, "transposon.fasta"), args.output)
            except Exception as e:
                print_error(f"Error processing {file_name}: {e}")
                continue
    elif os.path.isfile(args.input) and any(args.input.endswith(ext) for ext in valid_extensions):
        file_name = os.path.splitext(os.path.basename(args.input))[0]
        print_subheader(f"Processing single file: {file_name}")
        try:
            file_process_function(file_name, os.path.dirname(args.input), os.path.join(databases_dir, "Database.csv"),
                                os.path.join(databases_dir, "oric.fna"), os.path.join(databases_dir, "orit.fna"),
                                os.path.join(databases_dir, "plasmidfinder.fasta"), os.path.join(databases_dir, "transposon.fasta"), args.output)
        except Exception as e:
            print_error(f"Error processing {file_name}: {e}")
    else:
        print_error("Invalid path or file type. Please provide a valid directory or file.")

    # Clean up
    print_subheader("Cleaning up temporary files")
    shutil.rmtree('tmp_files', ignore_errors=True)
    shutil.rmtree('makedb_folder', ignore_errors=True)
    print_success("Temporary files removed")
    
    print_header("Plasmid Annotation Pipeline Completed")
    print_success(f"All output files are available in: {args.output}")

def cli():
    """Command line interface for the package."""
    main()

if __name__ == "__main__":
    main()