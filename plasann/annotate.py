import argparse
import os
from .essential_annotation import *
from .database_download import *
from .draw_plasmid import *
from . import blast_search as blast
import glob
import pandas as pd
import shutil
import sys
import platform
from Bio import SeqIO, Entrez
from Bio.Data import CodonTable
from collections import defaultdict
import uuid
import shutil
import time  
import threading
import uuid

# Version information
__version__ = "1.1.0"
__author__ = "Habibul Islam"
__email__ = "hislam2@ur.rochester.edu"

CDS_dataframe = pd.DataFrame()

def check_disk_space(required_gb=1):
    """Check available disk space"""
    try:
        free_bytes = shutil.disk_usage('.').free
        free_gb = free_bytes / (1024**3)
        
        if free_gb < required_gb:
            raise Exception(f"Insufficient disk space: {free_gb:.1f}GB available, {required_gb}GB required")
        return True
    except Exception as e:
        print(f"⚠️ Warning: Could not check disk space: {e}")
        return True  # Continue anyway

def print_version():
    """Print PlasAnn version information"""
    print(f"""
PlasAnn v{__version__}
Comprehensive Plasmid Annotation Pipeline

Author: {__author__}
Email: {__email__}
Python: {sys.version.split()[0]}
Platform: {platform.system()} {platform.release()}

For help: PlasAnn -h
For dependency check: PlasAnn --check-deps
""")

def print_dependency_status():
    """Print status of all external dependencies"""
    required_tools = {
        'makeblastdb': 'BLAST+',
        'blastn': 'BLAST+', 
        'blastx': 'BLAST+',
        'blastp': 'BLAST+',  # Added for UniProt BLAST
        'prodigal': 'Prodigal',
        'cmscan': 'Infernal',
        'cmpress': 'Infernal'
    }
    
    print(f"🔍 PlasAnn v{__version__} - External Dependency Status:")
    all_found = True
    
    for tool, package in required_tools.items():
        if shutil.which(tool):
            print(f"   ✅ {tool} ({package}) - Found")
        else:
            print(f"   ❌ {tool} ({package}) - Not found")
            all_found = False
    
    if all_found:
        print("\n✅ All external dependencies are installed!")
        print("🚀 PlasAnn is ready to use!")
    else:
        print("\n❌ Some dependencies are missing. Install them using:")
        system = platform.system().lower()
        if system == "darwin":
            print("   brew install blast prodigal infernal")
        elif system == "linux":
            print("   sudo apt install ncbi-blast+ prodigal infernal")
        else:
            print("   conda install -c bioconda blast prodigal infernal")
        print("\nThen run: PlasAnn --check-deps")

def check_external_tools():
    """Check if required external tools are installed"""
    required_tools = ['makeblastdb', 'blastn', 'blastx', 'blastp', 'prodigal', 'cmscan', 'cmpress']
    missing = [tool for tool in required_tools if not shutil.which(tool)]
    
    if missing:
        print("❌ Missing required external tools:")
        for tool in missing:
            print(f"   - {tool}")
        print("\nInstall with: conda install -c bioconda blast prodigal infernal")
        print("Or check: PlasAnn --check-deps")
        sys.exit(1)

def get_output_name(input_path, custom_name):
    if custom_name:
        return custom_name
    else:
        base = os.path.basename(input_path)
        return os.path.splitext(base)[0]

def validate_input(input_path, input_type):
    """Validate input file/folder based on specified type"""
    
    if not os.path.exists(input_path):
        print(f"❌ Error: Input path '{input_path}' does not exist.")
        sys.exit(1)
    
    fasta_extensions = ['.fasta', '.fa', '.fsa', '.fna']
    genbank_extensions = ['.gb', '.gbk', '.genbank']
    
    if input_type == "auto":
        if os.path.isfile(input_path):
            print("❌ Error: Auto detection requires a folder input, not a single file.")
            print("   For single files, please specify -t fasta or -t genbank explicitly.")
            sys.exit(1)
        
        elif os.path.isdir(input_path):
            fasta_files = []
            gb_files = []
            
            for ext in fasta_extensions:
                fasta_files.extend(glob.glob(os.path.join(input_path, f'*{ext}')))
            for ext in genbank_extensions:
                gb_files.extend(glob.glob(os.path.join(input_path, f'*{ext}')))
            
            total_files = len(fasta_files) + len(gb_files)
            
            if total_files == 0:
                print(f"❌ Error: Input folder '{input_path}' contains no FASTA or GenBank files.")
                print(f"   Expected FASTA extensions: {', '.join(fasta_extensions)}")
                print(f"   Expected GenBank extensions: {', '.join(genbank_extensions)}")
                sys.exit(1)
            else:
                print(f"✅ Auto-detection found:")
                if fasta_files:
                    print(f"   📄 {len(fasta_files)} FASTA file(s)")
                if gb_files:
                    print(f"   📄 {len(gb_files)} GenBank file(s)")
                print(f"   🎯 Total: {total_files} files to process")
    
    elif input_type == "fasta":
        if os.path.isfile(input_path):
            file_ext = os.path.splitext(input_path)[1].lower()
            if file_ext not in fasta_extensions:
                print(f"❌ Error: Input file '{input_path}' does not have a valid FASTA extension.")
                print(f"   Valid FASTA extensions: {', '.join(fasta_extensions)}")
                print(f"   Your file extension: {file_ext}")
                sys.exit(1)
        
        elif os.path.isdir(input_path):
            fasta_files = []
            for ext in fasta_extensions:
                fasta_files.extend(glob.glob(os.path.join(input_path, f'*{ext}')))
            
            if not fasta_files:
                print(f"❌ Error: Input folder '{input_path}' contains no FASTA files.")
                print(f"   Expected FASTA extensions: {', '.join(fasta_extensions)}")
                sys.exit(1)
            else:
                print(f"✅ Found {len(fasta_files)} FASTA file(s) in input folder.")
    
    elif input_type == "genbank":
        if os.path.isfile(input_path):
            file_ext = os.path.splitext(input_path)[1].lower()
            if file_ext not in genbank_extensions:
                print(f"❌ Error: Input file '{input_path}' does not have a valid GenBank extension.")
                print(f"   Valid GenBank extensions: {', '.join(genbank_extensions)}")
                print(f"   Your file extension: {file_ext}")
                sys.exit(1)
        
        elif os.path.isdir(input_path):
            gb_files = []
            for ext in genbank_extensions:
                gb_files.extend(glob.glob(os.path.join(input_path, f'*{ext}')))
            
            if not gb_files:
                print(f"❌ Error: Input folder '{input_path}' contains no GenBank files.")
                print(f"   Expected GenBank extensions: {', '.join(genbank_extensions)}")
                sys.exit(1)
            else:
                print(f"✅ Found {len(gb_files)} GenBank file(s) in input folder.")

def run_on_single_fasta(input_path, output_base, custom_name=None):
    """Process a single FASTA file"""
    global CDS_dataframe
    
    output_name = get_output_name(input_path, custom_name)
    output_folder = os.path.join(output_base, output_name)
    os.makedirs(output_folder, exist_ok=True)

    record = SeqIO.read(input_path, "fasta")
    fasta_length = len(record.seq)
    fasta_sequence = str(record.seq)

    df = CDS_fasta_using_prodigal(input_path, output_name) 
    CDS_dataframe = pd.concat([CDS_dataframe, df], ignore_index=True)
    #CDS_dataframe.to_csv('pyrodigal_output.csv', index=False)

    return CDS_dataframe, output_folder, output_name, fasta_length, fasta_sequence

def run_on_single_genbank(input_path, output_base, use_overwrite=False, custom_name=None):
    """Process a single GenBank file"""
    global CDS_dataframe
    extract_func = CDS_genbank_overwrite if use_overwrite else CDS_genbank_retain
    
    # Show processing mode
    if use_overwrite:
        print("🔄 Using Prodigal for gene prediction (--overwrite mode)")
    else:
        print("📋 Using GenBank annotations (--retain mode)")
    
    tmp_dir = "temp_dir"
    os.makedirs(tmp_dir, exist_ok=True)

    output_name = get_output_name(input_path, custom_name)
    output_folder = os.path.join(output_base, output_name)
    os.makedirs(output_folder, exist_ok=True)

    df, fasta_sequence, fasta_length = extract_func(input_path, output_name, tmp_dir)
    CDS_dataframe = pd.concat([CDS_dataframe, df], ignore_index=True)

    return CDS_dataframe, output_folder, output_name, fasta_length, fasta_sequence

'''def process_single_file_complete_pipeline(input_path, output_base, file_type, custom_name=None, 
                                         overwrite=False, shared_session=False, 
                                         uniprot_blast=False, uniprot_tsv=None, min_identity=30):
    """Process a single file through the complete annotation pipeline"""
    global CDS_dataframe
    
    # Reset global CDS_dataframe for each file
    CDS_dataframe = pd.DataFrame()
    
    print(f"🧬 Processing {file_type.upper()} file: {os.path.basename(input_path)}")
    
    # Use unique temp directories for each file when in shared session
    if shared_session:
        unique_id = str(uuid.uuid4())[:8]
        blast_temp_dir = f"temp_dir_blast_{unique_id}"
        temp_dir = f"temp_dir_{unique_id}"
    else:
        blast_temp_dir = "temp_dir_blast"
        temp_dir = "temp_dir"
    
    try:
        # Step 1: Extract CDS
        if file_type == "fasta":
            result = run_on_single_fasta(input_path, output_base, custom_name)
        elif file_type == "genbank":
            result = run_on_single_genbank(input_path, output_base, overwrite, custom_name)
        
        if result[0] is None or len(result[0]) == 0:
            print(f"❌ Failed to extract CDS from {input_path}")
            return False
            
        CDS_dataframe, output_folder, output_name, length, sequence = result
        
        # ✨ NEW: Check for overwrite mode without sequence
        if file_type == "genbank" and overwrite and sequence is None:
            print("❌ Add fasta file as input or use --retain option")
            return False
        
        # Step 2: Download and prepare databases (only if not shared session)
        if not shared_session:
            print("📥 Preparing databases...")
            download_database()
            prepare_blast_database()
        
        # Step 3: Run annotation pipeline
        print("🔍 Running BLAST annotation...")
        blast_db_prefix = "database_blast/translations_db"
        orit_fastas = ["Database/orit.fna", "Database/oriT_TNcentral.fasta"]
        rep_fasta = "Database/plasmidfinder.fasta"
        transposon_fastas = ["Database/tncentral_cleaned.fa","Database/transposon.fasta"]
        rfam_cm_path = "Database/Rfam.cm"
        
        blast.perform_blast_multiprocessing(CDS_dataframe, blast_db_prefix, blast_temp_dir)
        blast.annotate_blast_results(blast_temp_dir, "Database/Database.csv")
        dataframe_after_blast = generate_orf_annotation(CDS_dataframe, blast_temp_dir)
        
        # ✨ NEW: Check if sequence is available for mobile element detection
        if sequence is None:
            print("⚠️ Sequence content is not present. skipping oriC, oriT, ncRNA, transposon features")
            
            # Skip all sequence-dependent analyses
            dataframe_after_transposon = dataframe_after_blast
            ncrna_df = pd.DataFrame()  # Empty ncRNA dataframe
            
        else:
            # Full pipeline with sequence-dependent analyses
            print("🎯 Predicting origins and mobile elements...")
            dataframe_after_oric = blast.predict_oriC(sequence, dataframe_after_blast, "Database/oric.fna")
            dataframe_after_orit = blast.predict_oriT(sequence, dataframe_after_oric, orit_fastas)
            dataframe_after_replicon = blast.predict_replicons(sequence, dataframe_after_orit, rep_fasta)
            dataframe_after_transposon = blast.predict_transposons(sequence, dataframe_after_replicon, transposon_fastas)
            
            print("🧬 Detecting ncRNAs...")
            ncrna_df = blast.run_infernal_on_sequence(sequence, rfam_cm_path)
        
        print("📊 Finalizing annotations...")
        
        final_dataframe = combine_features(dataframe_after_transposon, ncrna_df, sequence)
        final_dataframe_fixed = Fixing_dataframe(final_dataframe, sequence)
        print_sequence_statistics(sequence, output_name, length)
        
        # ✨ Intergenic gene detection logic
        if sequence is not None:
            if file_type == "genbank" and not overwrite:
                # GenBank retain mode: trust original annotations, skip intergenic detection
                print("📋 Retaining original GenBank annotations - skipping intergenic gene detection")
            else:
                # FASTA files OR GenBank overwrite mode: run intergenic detection
                final_dataframe_fixed = blast.detect_intergenic_genes_simple(final_dataframe_fixed, sequence, blast_temp_dir)
        # If sequence is None, intergenic detection is skipped anyway
        
        # ✨ Optional UniProt BLAST annotation
        if uniprot_blast and uniprot_tsv:
            # Check if UniProt TSV exists after database download
            if os.path.exists(uniprot_tsv):
                print("🧬 Running optional UniProt BLAST annotation...")
                print("⚠️ This step may take several minutes...")
                final_dataframe_fixed = blast.uniprot_blast(final_dataframe_fixed, uniprot_tsv, min_identity)
                print("✅ UniProt BLAST annotation completed!")
            else:
                print(f"❌ UniProt TSV file not found: {uniprot_tsv}")
                print("   Make sure uniprot_plasmids.tsv is included in your database package")
                print("   Continuing without UniProt annotation...")
        elif uniprot_blast and not uniprot_tsv:
            print("❌ UniProt BLAST requested but no TSV file configured")
            print("   Continuing without UniProt annotation...")
        
        # Step 4: Write outputs
        output_csv = os.path.join(output_folder, f"{output_name}_annotations.csv")
        final_dataframe_fixed.to_csv(output_csv, index=False)

        
        
        # If UniProt was used, save a separate file for comparison
        if uniprot_blast and uniprot_tsv and os.path.exists(uniprot_tsv):
            output_csv_uniprot = os.path.join(output_folder, f"{output_name}_annotations_with_uniprot.csv")
            final_dataframe_fixed.to_csv(output_csv_uniprot, index=False)
            print(f"📊 UniProt-enhanced annotations saved: {output_csv_uniprot}")
        
        # ✨ Create GenBank file (with sequence handling)
        output_genbank = os.path.join(output_folder, f"{output_name}_genbank.gbk")
        description = "Annotated " + output_name
        write_genbank_from_annotation(final_dataframe_fixed, sequence, output_genbank, output_name, output_name, description)
        
        print("🎨 Creating plasmid map...")
        output_map = os.path.join(output_folder, f"{output_name}_map.png")
        draw_plasmid_map_from_genbank_file(output_genbank, output_map, output_name)
        
        # Clean up temp directories (but NOT shared database_blast in shared session)
        if shared_session:
            # Only clean up this file's unique temp dirs
            shutil.rmtree(blast_temp_dir, ignore_errors=True)
            shutil.rmtree(temp_dir, ignore_errors=True)
        else:
            # Clean up everything (single file mode)
            shutil.rmtree('database_blast', ignore_errors=True)
            shutil.rmtree('temp_dir_blast', ignore_errors=True)
            shutil.rmtree('temp_dir', ignore_errors=True)
            
        # Clean up UniProt BLAST database if it was created
        uniprot_db_dir = "uniprot_blast_db"
        if os.path.exists(uniprot_db_dir):
            shutil.rmtree(uniprot_db_dir, ignore_errors=True)
        
        print(f"✅ Completed processing {os.path.basename(input_path)}")
        print(f"   📁 Output folder: {output_folder}")
        
        # ✨ Show processing summary
        if sequence is None:
            cds_count = len(final_dataframe_fixed[final_dataframe_fixed.get('feature type', '') == 'CDS'])
            print(f"   📊 CDS-only mode: {cds_count} genes annotated")
            print(f"   ⚠️ Limited analysis due to missing genomic sequence")
        
        return True
        
    except Exception as e:
        print(f"❌ Error processing {input_path}: {e}")
        # Clean up on error
        if shared_session:
            shutil.rmtree(blast_temp_dir, ignore_errors=True)
            shutil.rmtree(temp_dir, ignore_errors=True)
        else:
            shutil.rmtree('database_blast', ignore_errors=True)
            shutil.rmtree('temp_dir_blast', ignore_errors=True)
            shutil.rmtree('temp_dir', ignore_errors=True)
        
        # Clean up UniProt BLAST database if it was created
        uniprot_db_dir = "uniprot_blast_db"
        if os.path.exists(uniprot_db_dir):
            shutil.rmtree(uniprot_db_dir, ignore_errors=True)
        return False'''



def process_single_file_complete_pipeline(input_path, output_base, file_type, custom_name=None, 
                                         overwrite=False, shared_session=False, 
                                         uniprot_blast=False, uniprot_tsv=None, min_identity=30):
    """Process a single file through the complete annotation pipeline - ENHANCED VERSION"""
    global CDS_dataframe
    
    # Reset global CDS_dataframe for each file
    CDS_dataframe = pd.DataFrame()
    
    print(f"🧬 Processing {file_type.upper()} file: {os.path.basename(input_path)}")
    
    # ✨ NEW: Pre-flight checks
    try:
        check_disk_space(required_gb=1)
    except Exception as e:
        print(f"❌ Pre-flight check failed: {e}")
        return False
    
    # Use unique temp directories for each file when in shared session
    if shared_session:
        unique_id = str(uuid.uuid4())[:8]
        blast_temp_dir = f"temp_dir_blast_{unique_id}"
        temp_dir = f"temp_dir_{unique_id}"
    else:
        blast_temp_dir = "temp_dir_blast"
        temp_dir = "temp_dir"
    
    try:
        # Step 1: Extract CDS with enhanced error handling
        if file_type == "fasta":
            result = run_on_single_fasta(input_path, output_base, custom_name)
        elif file_type == "genbank":
            result = run_on_single_genbank(input_path, output_base, overwrite, custom_name)
        
        # ✨ ENHANCED: Check for various failure modes
        if result is None:
            print(f"❌ Failed to process input file {input_path}")
            return False
        
        if len(result) != 5:
            print(f"❌ Unexpected result format from {input_path}")
            return False
            
        CDS_dataframe, output_folder, output_name, length, sequence = result
        
        # ✨ ENHANCED: Check CDS extraction results
        if CDS_dataframe is None:
            print(f"❌ CDS extraction returned None for {input_path}")
            return False
            
        if len(CDS_dataframe) == 0:
            print(f"⚠️ No genes detected in {input_path}")
            print("   This could mean:")
            print("   • Very short sequence (< 90bp)")
            print("   • Poor sequence quality")
            print("   • No open reading frames found")
            # Continue with empty dataframe - let user decide
        
        # ✨ NEW: Sequence quality checks
        if sequence is not None:
            # Check for high N content
            n_content = sequence.upper().count('N') / len(sequence) if len(sequence) > 0 else 0
            if n_content > 0.5:
                print(f"⚠️ High ambiguous nucleotide content ({n_content:.1%}) may affect results")
            
            # Check minimum length for meaningful analysis
            if len(sequence) < 200:
                print(f"⚠️ Very short sequence ({len(sequence)}bp) - results may be limited")
        
        # ✨ NEW: Check for overwrite mode without sequence
        if file_type == "genbank" and overwrite and sequence is None:
            print("❌ GenBank overwrite mode requires sequence data")
            print("   Solutions:")
            print("   • Use --retain mode instead of --overwrite")
            print("   • Provide FASTA file as input")
            print("   • Ensure GenBank file contains sequence data")
            return False
        
        # Step 2: Download and prepare databases (only if not shared session)
        if not shared_session:
            print("📥 Preparing databases...")
            try:
                download_database()
                prepare_blast_database()
            except Exception as e:
                print(f"❌ Database preparation failed: {e}")
                return False
        
        # ✨ ENHANCED: Validate database before BLAST
        blast_db_prefix = "database_blast/translations_db"
        if not os.path.exists(f"{blast_db_prefix}.pin"):
            print(f"❌ BLAST database not found: {blast_db_prefix}")
            return False
        
        # Step 3: Run annotation pipeline with enhanced error handling
        if len(CDS_dataframe) > 0:
            print("🔍 Running BLAST annotation...")
            try:
                blast.perform_blast_multiprocessing(CDS_dataframe, blast_db_prefix, blast_temp_dir)
                blast.annotate_blast_results(blast_temp_dir, "Database/Database.csv")
                dataframe_after_blast = generate_orf_annotation(CDS_dataframe, blast_temp_dir)
            except Exception as e:
                print(f"❌ BLAST annotation failed: {e}")
                print("   Continuing with basic gene predictions only...")
                # Create basic annotations without BLAST results
                dataframe_after_blast = CDS_dataframe.copy()
                for col in ["Gene Name", "Product", "Category"]:
                    if col not in dataframe_after_blast.columns:
                        dataframe_after_blast[col] = "ORF" if col == "Gene Name" else "Open reading frame"
        else:
            print("⏭️ No genes found, skipping BLAST annotation")
            dataframe_after_blast = pd.DataFrame()
        
        # ✨ ENHANCED: Check if sequence is available for mobile element detection
        if sequence is None or len(sequence) == 0:
            print("⚠️ Sequence content not available - skipping sequence-dependent analyses:")
            print("   • Origin of replication (oriC) prediction")
            print("   • Origin of transfer (oriT) detection") 
            print("   • Transposon and mobile element detection")
            print("   • ncRNA detection")
            print("   • Intergenic gene discovery")
            
            dataframe_after_transposon = dataframe_after_blast
            ncrna_df = pd.DataFrame()
            
        else:
            # Full pipeline with sequence-dependent analyses
            print("🎯 Predicting origins and mobile elements...")
            try:
                orit_fastas = ["Database/orit.fna", "Database/oriT_TNcentral.fasta"]
                rep_fasta = "Database/plasmidfinder.fasta"
                transposon_fastas = ["Database/tncentral_cleaned.fa","Database/transposon.fasta"]
                rfam_cm_path = "Database/Rfam.cm"
                
                dataframe_after_oric = blast.predict_oriC(sequence, dataframe_after_blast, "Database/oric.fna")
                dataframe_after_orit = blast.predict_oriT(sequence, dataframe_after_oric, orit_fastas)
                dataframe_after_replicon = blast.predict_replicons(sequence, dataframe_after_orit, rep_fasta)
                dataframe_after_transposon = blast.predict_transposons(sequence, dataframe_after_replicon, transposon_fastas)
                
                print("🧬 Detecting ncRNAs...")
                ncrna_df = blast.run_infernal_on_sequence(sequence, rfam_cm_path)
                
            except Exception as e:
                print(f"⚠️ Mobile element detection failed: {e}")
                print("   Continuing with gene annotations only...")
                dataframe_after_transposon = dataframe_after_blast
                ncrna_df = pd.DataFrame()
        
        print("📊 Finalizing annotations...")
        
        final_dataframe = combine_features(dataframe_after_transposon, ncrna_df, sequence)
        final_dataframe_fixed = Fixing_dataframe(final_dataframe, sequence)
        print_sequence_statistics(sequence, output_name, length)
        
        # ✨ ENHANCED: Intergenic gene detection with better error handling
        if sequence is not None and len(sequence) > 0:
            if file_type == "genbank" and not overwrite:
                print("📋 Retaining original GenBank annotations - skipping intergenic gene detection")
            else:
                try:
                    print("🔍 Analyzing intergenic regions...")
                    final_dataframe_fixed = blast.detect_intergenic_genes_simple(final_dataframe_fixed, sequence, blast_temp_dir)
                except Exception as e:
                    print(f"⚠️ Intergenic gene detection failed: {e}")
                    print("   Continuing with current annotations...")
        
        # ✨ Enhanced UniProt BLAST with better error handling
        if uniprot_blast and uniprot_tsv:
            if os.path.exists(uniprot_tsv):
                try:
                    print("🧬 Running optional UniProt BLAST annotation...")
                    print("⚠️ This step may take several minutes...")
                    final_dataframe_fixed = blast.uniprot_blast(final_dataframe_fixed, uniprot_tsv, min_identity)
                    print("✅ UniProt BLAST annotation completed!")
                except Exception as e:
                    print(f"⚠️ UniProt BLAST failed: {e}")
                    print("   Continuing without UniProt annotations...")
            else:
                print(f"❌ UniProt TSV file not found: {uniprot_tsv}")
                print("   Continuing without UniProt annotation...")
        elif uniprot_blast and not uniprot_tsv:
            print("❌ UniProt BLAST requested but no TSV file configured")
        
        # Step 4: Write outputs with error handling
        try:
            output_csv = os.path.join(output_folder, f"{output_name}_annotations.csv")
            final_dataframe_fixed.to_csv(output_csv, index=False)
            
            if uniprot_blast and uniprot_tsv and os.path.exists(uniprot_tsv):
                output_csv_uniprot = os.path.join(output_folder, f"{output_name}_annotations_with_uniprot.csv")
                final_dataframe_fixed.to_csv(output_csv_uniprot, index=False)
                print(f"📊 UniProt-enhanced annotations saved: {output_csv_uniprot}")
            
        except Exception as e:
            print(f"❌ Failed to write CSV outputs: {e}")
            return False
        
        # ✨ ENHANCED: GenBank and map generation with better error handling
        try:
            output_genbank = os.path.join(output_folder, f"{output_name}_genbank.gbk")
            description = "Annotated " + output_name
            write_genbank_from_annotation(final_dataframe_fixed, sequence, output_genbank, output_name, output_name, description)
            
            print("🎨 Creating plasmid map...")
            output_map = os.path.join(output_folder, f"{output_name}_map.png")
            
            map_success = draw_plasmid_map_from_genbank_file(output_genbank, output_map, output_name)
            if not map_success:
                print("⚠️ Plasmid map generation failed, but analysis completed successfully")
                
        except Exception as e:
            print(f"⚠️ GenBank/map generation failed: {e}")
            print("   CSV annotations are still available")
        
        # Clean up temp directories (but NOT shared database_blast in shared session)
        if shared_session:
            shutil.rmtree(blast_temp_dir, ignore_errors=True)
            shutil.rmtree(temp_dir, ignore_errors=True)
        else:
            shutil.rmtree('database_blast', ignore_errors=True)
            shutil.rmtree('temp_dir_blast', ignore_errors=True)
            shutil.rmtree('temp_dir', ignore_errors=True)
            
        # Clean up UniProt BLAST database if it was created
        uniprot_db_dir = "uniprot_blast_db"
        if os.path.exists(uniprot_db_dir):
            shutil.rmtree(uniprot_db_dir, ignore_errors=True)
        
        print(f"✅ Completed processing {os.path.basename(input_path)}")
        print(f"   📁 Output folder: {output_folder}")
        
        # ✨ ENHANCED: Processing summary
        if sequence is None:
            cds_count = len(final_dataframe_fixed[final_dataframe_fixed.get('feature type', '') == 'CDS'])
            print(f"   📊 CDS-only mode: {cds_count} genes annotated")
            print(f"   ⚠️ Limited analysis due to missing genomic sequence")
        else:
            total_features = len(final_dataframe_fixed)
            cds_count = len(final_dataframe_fixed[final_dataframe_fixed.get('feature type', '') == 'CDS'])
            print(f"   📊 Complete analysis: {total_features} total features ({cds_count} genes)")
        
        return True
        
    except Exception as e:
        print(f"❌ Error processing {input_path}: {e}")
        import traceback
        print("   Full error details:")
        traceback.print_exc()
        
        # Enhanced cleanup on error
        cleanup_dirs = [blast_temp_dir, temp_dir, 'database_blast', 'temp_dir_blast', 'uniprot_blast_db']
        for cleanup_dir in cleanup_dirs:
            if os.path.exists(cleanup_dir):
                shutil.rmtree(cleanup_dir, ignore_errors=True)
        
        return False

def process_folder(input_folder, output_base, file_type, custom_name=None, overwrite=False,
                  uniprot_blast=False, uniprot_tsv=None, min_identity=30):
    """Process all files in a folder of specified type"""
    
    if file_type == "fasta":
        extensions = ['.fasta', '.fa', '.fsa', '.fna']
    elif file_type == "genbank":
        extensions = ['.gb', '.gbk', '.genbank']
    
    # Find all files of the specified type
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(input_folder, f'*{ext}')))
    
    total_files = len(files)
    successful_files = 0
    
    print(f"\n🎯 Processing {total_files} {file_type.upper()} files from folder...")
    
    # Show UniProt status
    if uniprot_blast:
        print(f"🧬 UniProt BLAST: Enabled (TSV: {os.path.basename(uniprot_tsv) if uniprot_tsv else 'default'})")
        print(f"   ⚠️ This will significantly increase processing time!")
    
    # Prepare databases ONCE for all files
    print("📥 Preparing shared databases for all files...")
    download_database()
    prepare_blast_database()
    
    # Process each file through pipeline (with shared databases)
    for i, file_path in enumerate(files, 1):
        print(f"\n{'='*60}")
        print(f"📁 File {i}/{total_files}: {os.path.basename(file_path)}")
        print(f"{'='*60}")
        
        # shared_session=True means don't rebuild databases, use unique temp dirs
        success = process_single_file_complete_pipeline(
            file_path, output_base, file_type, custom_name, overwrite, 
            shared_session=True, uniprot_blast=uniprot_blast, 
            uniprot_tsv=uniprot_tsv, min_identity=min_identity
        )
        if success:
            successful_files += 1
    
    # Clean up shared databases after all files are done
    print("🧹 Cleaning up shared databases...")
    shutil.rmtree('database_blast', ignore_errors=True)
    
    # Clean up UniProt BLAST database if it was used
    uniprot_db_dir = "uniprot_blast_db"
    if os.path.exists(uniprot_db_dir):
        shutil.rmtree(uniprot_db_dir, ignore_errors=True)
    
    print(f"\n🎉 Folder processing complete!")
    print(f"   ✅ Successfully processed: {successful_files}/{total_files} files")
    
    if uniprot_blast and successful_files > 0:
        print(f"   🧬 UniProt BLAST annotations were included in results")
    
    if successful_files == 0:
        print("❌ No files were successfully processed.")
        sys.exit(1)


def process_auto_folder(input_folder, output_base, custom_name=None, overwrite=False,
                       uniprot_blast=False, uniprot_tsv=None, min_identity=30):
    """Process folder containing mixed FASTA and GenBank files"""
    
    fasta_extensions = ['.fasta', '.fa', '.fsa', '.fna']
    genbank_extensions = ['.gb', '.gbk', '.genbank']
    
    # Find all files
    fasta_files = []
    gb_files = []
    
    for ext in fasta_extensions:
        fasta_files.extend(glob.glob(os.path.join(input_folder, f'*{ext}')))
    for ext in genbank_extensions:
        gb_files.extend(glob.glob(os.path.join(input_folder, f'*{ext}')))
    
    total_files = len(fasta_files) + len(gb_files)
    successful_files = 0
    
    print(f"\n🎯 Starting auto-detection pipeline for {total_files} files...")
    
    # Show GenBank processing method if there are GenBank files
    if gb_files:
        method = "Prodigal (--overwrite)" if overwrite else "GenBank annotations (--retain)"
        print(f"   📄 GenBank processing method: {method}")
    
    # Show UniProt status
    if uniprot_blast:
        print(f"   🧬 UniProt BLAST: Enabled (TSV: {os.path.basename(uniprot_tsv) if uniprot_tsv else 'default'})")
        print(f"   ⚠️ This will significantly increase processing time!")
    
    # Prepare databases ONCE for all files
    print("📥 Preparing shared databases for all files...")
    download_database()
    prepare_blast_database()
    
    # Process FASTA files (with shared databases)
    for i, fasta_file in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"📁 FASTA file {i}/{len(fasta_files)}: {os.path.basename(fasta_file)}")
        print(f"{'='*60}")
        
        success = process_single_file_complete_pipeline(
            fasta_file, output_base, "fasta", custom_name, False, 
            shared_session=True, uniprot_blast=uniprot_blast, 
            uniprot_tsv=uniprot_tsv, min_identity=min_identity
        )
        if success:
            successful_files += 1
    
    # Process GenBank files with the specified method (with shared databases)
    for i, gb_file in enumerate(gb_files, 1):
        print(f"\n{'='*60}")
        print(f"📁 GenBank file {i}/{len(gb_files)}: {os.path.basename(gb_file)}")
        print(f"{'='*60}")
        
        success = process_single_file_complete_pipeline(
            gb_file, output_base, "genbank", custom_name, overwrite, 
            shared_session=True, uniprot_blast=uniprot_blast, 
            uniprot_tsv=uniprot_tsv, min_identity=min_identity
        )
        if success:
            successful_files += 1
    
    # Clean up shared databases after all files are done
    print("🧹 Cleaning up shared databases...")
    shutil.rmtree('database_blast', ignore_errors=True)
    
    # Clean up UniProt BLAST database if it was used
    uniprot_db_dir = "uniprot_blast_db"
    if os.path.exists(uniprot_db_dir):
        shutil.rmtree(uniprot_db_dir, ignore_errors=True)
    
    print(f"\n🎉 Auto-detection complete!")
    print(f"   ✅ Successfully processed: {successful_files}/{total_files} files")
    
    if uniprot_blast and successful_files > 0:
        print(f"   🧬 UniProt BLAST annotations were included in results")
    
    if successful_files == 0:
        print("❌ No files were successfully processed.")
        sys.exit(1)


def CDS_genbank_retain(genbank_file, output_name, email="hislam2@ur.rochester.edu"):
    """Extract CDS using GenBank annotations and fallback to reverse translation if sequence is missing."""
    Entrez.email = email
    cds_info = []
    records = list(SeqIO.parse(genbank_file, "genbank"))
    
    full_seq = None
    accession = None

    if records:
        accession = records[0].id

    # Try fetching the sequence from NCBI
    try:
        print(f"🔗 Fetching full FASTA for accession: {accession}")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        fasta_record = SeqIO.read(handle, "fasta")
        handle.close()
        full_seq = fasta_record.seq.lower()
    except Exception as e:
        print(f"⚠️ Could not fetch sequence for {accession}: Using Translation sequence from GenBank input")
        full_seq = None

    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue

            if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
                start = int(feature.location.parts[0].start) + 1
                end = int(feature.location.parts[-1].end)
            else:
                start = int(feature.location.start) + 1
                end = int(feature.location.end)

            strand = '+' if feature.location.strand == 1 else '-'
            seq = None

            if full_seq:
                try:
                    seq = feature.extract(full_seq)
                    if feature.location.strand == -1:
                        seq = seq.reverse_complement()
                    seq = str(seq).lower()
                except Exception as e:
                    print(f"⚠️ Failed to extract sequence for {start}..{end}: {e}")
                    seq = None

            # Fallback: reverse translate from protein
            if not seq and "translation" in feature.qualifiers:
                protein_seq = feature.qualifiers["translation"][0]
                bacterial_table = CodonTable.unambiguous_dna_by_id[11]
                aa_to_codon = defaultdict(list)
                for codon, aa in bacterial_table.forward_table.items():
                    aa_to_codon[aa].append(codon)
                codon_map = {aa: sorted(codons)[0] for aa, codons in aa_to_codon.items()}
                seq = ""
                for aa in protein_seq:
                    if aa == '*':
                        break
                    seq += codon_map.get(aa, 'NNN')

            if seq:
                cds_info.append({
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "Length": len(seq),
                    "Sequence": seq
                })

    if cds_info:
        df = pd.DataFrame(cds_info)
        print(f"✅ Extracted {len(df)} CDS entries from {output_name} plasmid")
    else:
        print("⚠️ No CDS entries extracted.")
        df = pd.DataFrame()

    return df, str(full_seq) if full_seq else None, len(full_seq) if full_seq else None

def main():
    parser = argparse.ArgumentParser(
        description=f"""
PlasAnn v{__version__} - Comprehensive Plasmid Annotation Pipeline

🧬 Features:
  • Gene prediction (Prodigal) and functional annotation (BLAST)
  • Mobile element detection (oriC, oriT, transposons, replicons)
  • ncRNA detection (Infernal/Rfam) and intergenic gene discovery
  • Optional UniProt BLAST enhancement for comprehensive annotation
  • Beautiful circular plasmid visualizations
  • Batch processing with auto-detection of mixed file types

📖 For detailed documentation, visit: https://github.com/your-repo/plasann
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📚 Examples:

  Basic Usage:
    PlasAnn -i plasmid.fasta -o results -t fasta
    PlasAnn -i plasmid.gb -o results -t genbank
    PlasAnn -i mixed_folder/ -o results -t auto
  
  Enhanced Annotation:
    PlasAnn -i plasmid.fasta -o results -t fasta --uniprot-blast
  
  GenBank Processing Modes:
    PlasAnn -i plasmid.gb -o results -t genbank --retain   # Use original annotations
    PlasAnn -i plasmid.gb -o results -t genbank --overwrite # Re-annotate with Prodigal
  
  Batch Processing:
    PlasAnn -i fasta_folder/ -o results -t fasta
    PlasAnn -i mixed_folder/ -o results -t auto --uniprot-blast
  
  System Check:
    PlasAnn --check-deps
    PlasAnn --version

🔧 Dependencies: BLAST+, Prodigal, Infernal
   Install with: conda install -c bioconda blast prodigal infernal

💡 Tip: Use --uniprot-blast for comprehensive protein annotation (slower but thorough)
        """
    )
    
    # Version and help arguments
    parser.add_argument("--version", action="store_true", 
                       help="Show PlasAnn version and exit")
    parser.add_argument("--check-deps", action="store_true", 
                       help="Check external dependency status and exit")
    
    # Required arguments
    parser.add_argument("-i", "--input", help="Input FASTA or GenBank file/folder")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-t", "--type", choices=["fasta", "genbank", "auto"], 
                       help="Input type: fasta, genbank, or auto (auto-detect from folder)")
    
    # Optional arguments
    parser.add_argument("-n", "--name", help="Custom name for output subfolder")
    parser.add_argument("--overwrite", action="store_true", 
                       help="Use Prodigal on GenBank sequence (ignore existing annotations)")
    parser.add_argument("--retain", action="store_true", 
                       help="Use GenBank annotations with fallback to translation (default)")
    
    # UniProt BLAST arguments
    parser.add_argument("--uniprot-blast", action="store_true", 
                       help="Run optional UniProt BLAST annotation (slow but comprehensive)")
    parser.add_argument("--uniprot-tsv", type=str, 
                       help="Path to UniProt TSV file (default: Database/uniprot_plasmids.tsv)")
    parser.add_argument("--min-identity", type=float, default=30,
                       help="Minimum identity percentage for UniProt BLAST hits (default: 30%%)")

    args = parser.parse_args()

    # Handle version check first
    if args.version:
        print_version()
        return

    # Handle dependency check
    if args.check_deps:
        print_dependency_status()
        return

    # Validate required arguments when NOT using --version or --check-deps
    if not args.input or not args.output or not args.type:
        parser.error("Arguments -i/--input, -o/--output, and -t/--type are required")

    # Handle UniProt TSV file path
    if args.uniprot_blast:
        # Use provided path or default to Database location
        if args.uniprot_tsv:
            uniprot_tsv_path = args.uniprot_tsv
        else:
            uniprot_tsv_path = "Database/uniprot_plasmids.tsv"
        
        # Check if file exists (will be available after database download)
        args.uniprot_tsv = uniprot_tsv_path
    else:
        args.uniprot_tsv = None

    # Validate input before proceeding
    validate_input(args.input, args.type)

    # Check dependencies before running
    check_external_tools()
    
    # Show configuration
    print(f"\n🔧 PlasAnn v{__version__} Configuration:")
    print(f"   📄 Input: {args.input}")
    print(f"   📁 Output: {args.output}")
    print(f"   🎯 Type: {args.type}")
    if args.uniprot_blast:
        print(f"   🧬 UniProt BLAST: Enabled")
        print(f"   📊 UniProt TSV: {args.uniprot_tsv}")
        print(f"   🎚️ Min Identity: {args.min_identity}%")
        print(f"   ⚠️  Note: UniProt BLAST will significantly increase runtime")
        # Check if file exists, but don't fail here (it might be downloaded with databases)
        if not os.path.exists(args.uniprot_tsv):
            print(f"   📥 UniProt TSV will be available after database download")
    else:
        print(f"   🧬 UniProt BLAST: Disabled")
    
    # Handle different input types
    if args.type == "auto":
        # Auto mode: determine GenBank processing method
        use_overwrite = args.overwrite and not args.retain
        process_auto_folder(args.input, args.output, args.name, use_overwrite,
                           args.uniprot_blast, args.uniprot_tsv, args.min_identity)
        
    elif os.path.isdir(args.input):
        # Folder mode: process all files of specified type
        if args.type == "genbank":
            use_overwrite = args.overwrite and not args.retain
            process_folder(args.input, args.output, args.type, args.name, use_overwrite,
                          args.uniprot_blast, args.uniprot_tsv, args.min_identity)
        else:
            # FASTA folder - overwrite flag doesn't apply
            process_folder(args.input, args.output, args.type, args.name, False,
                          args.uniprot_blast, args.uniprot_tsv, args.min_identity)
        
    else:
        # Single file mode: process one file (shared_session=False means build databases)
        if args.type == "genbank":
            use_overwrite = args.overwrite and not args.retain
            success = process_single_file_complete_pipeline(
                args.input, args.output, args.type, args.name, use_overwrite, 
                shared_session=False, uniprot_blast=args.uniprot_blast, 
                uniprot_tsv=args.uniprot_tsv, min_identity=args.min_identity
            )
        else:
            # FASTA file - overwrite flag doesn't apply
            success = process_single_file_complete_pipeline(
                args.input, args.output, args.type, args.name, False, 
                shared_session=False, uniprot_blast=args.uniprot_blast, 
                uniprot_tsv=args.uniprot_tsv, min_identity=args.min_identity
            )
            
        if not success:
            print("❌ Processing failed")
            sys.exit(1)

def cli_main():
    """Entry point for pip-installed PlasAnn command"""
    try:
        main()
    except KeyboardInterrupt:
        print("\n❌ Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)
    
if __name__ == "__main__":
    main()