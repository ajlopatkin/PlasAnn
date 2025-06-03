import os
import pandas as pd
from pathlib import Path
import shutil
from Bio import SeqIO
import re
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import tempfile
from Bio.Seq import Seq
import hashlib
import subprocess
import multiprocessing
from pathlib import Path


def validate_sequence_for_blast(seq, min_length=90):
    """Validate sequence before BLAST"""
    if pd.isna(seq) or not seq.strip():
        return False, "Empty sequence"
    
    seq_clean = seq.strip().upper()
    
    # Check minimum length
    if len(seq_clean) < min_length:
        return False, f"Too short ({len(seq_clean)}bp < {min_length}bp)"
    
    # Check for valid nucleotides
    valid_bases = set('ACGTRYSWKMBDHVN')
    invalid_bases = set(seq_clean) - valid_bases
    if invalid_bases:
        return False, f"Invalid nucleotides: {invalid_bases}"
    
    # Check for too many N's
    n_content = seq_clean.count('N') / len(seq_clean)
    if n_content > 0.8:
        return False, f"Too many ambiguous bases ({n_content:.1%})"
    
    return True, "Valid"

def run_blast_for_row(idx, seq, tmp_dir, blast_db_prefix, blast_header):
    """Enhanced BLAST execution with better error handling"""
    
    # ✨ NEW: Validate sequence first
    is_valid, reason = validate_sequence_for_blast(seq)
    if not is_valid:
        print(f"⚠️ Skipping ORF {idx}: {reason}")
        return

    tmp_dir = Path(tmp_dir)
    query_file = tmp_dir / f"query_{idx}.fasta"
    result_file = tmp_dir / f"blast_result_{idx}.csv"

    try:
        # Write query FASTA
        with open(query_file, "w") as f:
            f.write(f">query_{idx}\n{seq.strip()}\n")

        # ✨ ENHANCED: Use subprocess with proper error handling
        blast_cmd = [
            "blastx", 
            "-query", str(query_file),
            "-db", blast_db_prefix,
            "-out", str(result_file),
            "-outfmt", "10",
            "-max_target_seqs", "10", 
            "-evalue", "1e-5",
            "-num_threads", "1"  # Single thread per process
        ]
        
        # Run BLAST with timeout
        result = subprocess.run(
            blast_cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout per sequence
            check=True
        )
        
        # ✨ NEW: Check if BLAST produced output
        if result_file.exists():
            with open(result_file, "r") as f:
                content = f.read().strip()
            
            # Distinguish between no hits vs failed BLAST
            if content:  # Has BLAST hits
                with open(result_file, "w") as f:
                    f.write(blast_header + "\n" + content)
            else:  # No hits found (but BLAST ran successfully)
                with open(result_file, "w") as f:
                    f.write(blast_header + "\n")  # Header only
        else:
            # BLAST failed to create output file
            print(f"⚠️ BLAST failed to create output for ORF {idx}")
            # Create empty file with header
            with open(result_file, "w") as f:
                f.write(blast_header + "\n")
                
    except subprocess.TimeoutExpired:
        print(f"⏰ BLAST timeout for ORF {idx} (>5 minutes)")
        # Create empty result file
        with open(result_file, "w") as f:
            f.write(blast_header + "\n")
            
    except subprocess.CalledProcessError as e:
        print(f"❌ BLAST failed for ORF {idx}: {e.stderr}")
        # Create empty result file  
        with open(result_file, "w") as f:
            f.write(blast_header + "\n")
            
    except Exception as e:
        print(f"❌ Unexpected error for ORF {idx}: {e}")
        # Create empty result file
        with open(result_file, "w") as f:
            f.write(blast_header + "\n")
    
    finally:
        # Clean up query file
        if query_file.exists():
            query_file.unlink()

def perform_blast_multiprocessing(CDS_dataframe, blast_db_prefix, tmp_dir, max_workers=8):
    """Enhanced multiprocessing BLAST with validation"""
    
    # ✨ NEW: Pre-flight checks
    if not os.path.exists(f"{blast_db_prefix}.pin"):
        raise Exception(f"BLAST database not found: {blast_db_prefix}")
    
    if len(CDS_dataframe) == 0:
        print("⚠️ No sequences to BLAST")
        return
    
    # ✨ NEW: Filter valid sequences
    valid_sequences = []
    for idx, row in CDS_dataframe.iterrows():
        seq = row.get("Sequence", "")
        is_valid, reason = validate_sequence_for_blast(seq)
        if is_valid:
            valid_sequences.append((idx, seq))
        else:
            print(f"⚠️ Skipping ORF {idx}: {reason}")
    
    if not valid_sequences:
        print("❌ No valid sequences for BLAST")
        return
    
    print(f"🔍 Running BLAST on {len(valid_sequences)}/{len(CDS_dataframe)} sequences")
    
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    blast_header = ",".join([
        "query_id", "subject_id", "identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end",
        "s_start", "s_end", "evalue", "bit_score"
    ])

    # ✨ ENHANCED: Adaptive worker count based on system and sequence count
    import multiprocessing
    system_cores = multiprocessing.cpu_count()
    adaptive_workers = min(max_workers, system_cores, len(valid_sequences))
    
    print(f"   Using {adaptive_workers} parallel BLAST processes")

    tasks = []
    completed = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=adaptive_workers) as executor:
        # Submit all tasks
        for idx, seq in valid_sequences:
            tasks.append(executor.submit(
                run_blast_for_row, idx, seq, str(tmp_dir), blast_db_prefix, blast_header
            ))

        # Process completed tasks with progress tracking
        for i, future in enumerate(as_completed(tasks)):
            try:
                future.result()
                completed += 1
            except Exception as e:
                print(f"❌ Task {i} failed: {e}")
                failed += 1
            
            # Progress indicator
            total_done = completed + failed
            if total_done % 10 == 0 or total_done == len(tasks):
                print(f"   Progress: {total_done}/{len(tasks)} sequences processed")

    print(f"✅ BLAST completed: {completed} successful, {failed} failed")
    
    # ✨ NEW: Validate output files
    expected_files = len(valid_sequences)
    actual_files = len(list(tmp_dir.glob("blast_result_*.csv")))
    
    if actual_files < expected_files:
        print(f"⚠️ Warning: Expected {expected_files} result files, found {actual_files}")

def annotate_blast_results(blast_results_dir, database_csv_path):
    """Enhanced BLAST result annotation with better error handling"""
    blast_results_dir = Path(blast_results_dir)
    
    # ✨ NEW: Validate database file
    if not os.path.exists(database_csv_path):
        raise Exception(f"Database CSV not found: {database_csv_path}")
    
    try:
        db = pd.read_csv(database_csv_path)
        if len(db) == 0:
            raise Exception("Database CSV is empty")
    except Exception as e:
        raise Exception(f"Could not read database CSV: {e}")

    columns_to_add = ["Gene Name", "Category", "Product", "COG_DESCRIPTION", "KEGG_Function"]
    
    # ✨ NEW: Check if required columns exist
    missing_cols = [col for col in columns_to_add if col not in db.columns]
    if missing_cols:
        print(f"⚠️ Missing database columns: {missing_cols}")
        # Add empty columns for missing ones
        for col in missing_cols:
            db[col] = ""

    blast_files = list(blast_results_dir.glob("blast_result_*.csv"))
    if not blast_files:
        print("⚠️ No BLAST result files found in:", blast_results_dir)
        return

    processed = 0
    failed = 0

    for blast_file in blast_files:
        try:
            # ✨ NEW: Check if file has content beyond header
            with open(blast_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) <= 1:  # Only header or empty
                continue  # Skip annotation for files with no hits
            
            blast_df = pd.read_csv(blast_file)
            
            if blast_df.empty:
                continue
            
            # ✨ NEW: Validate subject_id column
            if 'subject_id' not in blast_df.columns:
                print(f"⚠️ No subject_id column in {blast_file.name}")
                continue
            
            # Convert to int and handle invalid values
            blast_df = blast_df.dropna(subset=['subject_id'])
            blast_df['subject_id'] = blast_df['subject_id'].astype(int)
            
            # ✨ NEW: Check for valid subject IDs
            valid_ids = blast_df['subject_id'].isin(db.index)
            if not valid_ids.any():
                print(f"⚠️ No valid subject IDs in {blast_file.name}")
                continue
            
            blast_df = blast_df[valid_ids]

            metadata = db.loc[blast_df["subject_id"], columns_to_add].reset_index(drop=True)
            annotated = pd.concat([blast_df, metadata], axis=1)
            annotated.to_csv(blast_file, index=False)
            processed += 1

        except Exception as e:
            print(f"❌ Failed to annotate {blast_file.name}: {e}")
            failed += 1

    print(f"✅ Annotated {processed} BLAST result files ({failed} failed)")

'''def predict_oriC(sequence,annotation_df,doric_fasta,min_length=200,location_buffer=50,window_size=500,sequence_id="plasmid_query"):

    seq = sequence.lower()
    disparity_curves = {
        "GC": {"positive": "g", "negative": "c"},
        "AT": {"positive": "a", "negative": "t"},
        "MK": {"positive": "a", "negative": "c"},
        "RY": {"positive": "a", "negative": "g"},
    }

    results = {}
    for name, groups in disparity_curves.items():
        curve = []
        pos = groups["positive"]
        neg = groups["negative"]
        for i in range(0, len(seq) - window_size + 1, window_size):
            window = seq[i:i + window_size]
            disparity = window.count(pos) - window.count(neg)
            curve.append(disparity)
        results[name] = np.cumsum(curve)

    normalized_curves = {k: (v - np.mean(v)) / np.std(v) for k, v in results.items()}
    combined_score = sum(normalized_curves.values())
    min_index = np.argmin(combined_score)
    center = (min_index * window_size + (min_index + 1) * window_size) // 2
    extract_start = max(0, center - 1000)
    extract_end = min(len(seq), center + 1000)
    oriC_seq = seq[extract_start:extract_end]

    dnaa_motifs = [ "ttatccaca", "ttatccaaa", "ttgtccaca", "ttatgcaca", "ttttccaca", "ttatccata", "ttatccacc"]
    motif_found = any(re.search(m, oriC_seq) for m in dnaa_motifs)

    # === Create temporary FASTA file ===
    temp_fasta = f"{sequence_id}_temp.fasta"
    with open(temp_fasta, "w") as f:
        f.write(f">{sequence_id}\n{seq}\n")

    # === Set up BLAST database folder ===
    db_folder = Path("oric_db_folder")
    db_folder.mkdir(exist_ok=True)
    doric_db_name = str(db_folder / "doric_temp_db")

    # === Run BLAST ===
    blast_output = "oric_blast_hits.csv"
    os.system(f"makeblastdb -in {doric_fasta} -dbtype nucl -out {doric_db_name}")
    os.system(
        f"blastn -query {temp_fasta} -db {doric_db_name} "
        f"-evalue 1e-5 -outfmt '10 qseqid sseqid pident length mismatch gapopen "
        f"qstart qend sstart send evalue bitscore stitle' "
        f"-max_target_seqs 10 -out {blast_output}"
    )

    # === Parse BLAST results ===
    header = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    doric_df = pd.DataFrame()
    if os.path.exists(blast_output):
        try:
            doric_df = pd.read_csv(blast_output, header=None, names=header)
            doric_df = doric_df[doric_df["length"] >= min_length].sort_values("qstart")
            filtered_hits = []
            last_end = -float("inf")
            for _, row in doric_df.iterrows():
                if row.qstart - last_end > location_buffer:
                    filtered_hits.append(row)
                    last_end = row.qend
            doric_df = pd.DataFrame(filtered_hits)
        except Exception:
            doric_df = pd.DataFrame()

    # === Append oriC entry ===
    oriC_entry = {
        "Gene Name": "oriV",
        "Product": "Predicted origin of replication",
        "Start": extract_start,
        "End": extract_end,
        "Strand": "+",
        "Category": "Origin of Replication",
        "Translation": oriC_seq.upper()
    }
    for col in annotation_df.columns:
        if col not in oriC_entry:
            oriC_entry[col] = ""

    annotation_df = pd.concat([annotation_df, pd.DataFrame([oriC_entry])], ignore_index=True)

    # === Clean up ===
    os.remove(temp_fasta)
    if os.path.exists(blast_output):
        os.remove(blast_output)
    if db_folder.exists():
        shutil.rmtree(db_folder)

    print("[INFO] oriC prediction complete.")
    return annotation_df'''

def predict_oriC(sequence, annotation_df, doric_fasta, min_length=100, location_buffer=50, 
                 window_size=500, sequence_id="plasmid_query", max_search_region=3000):
    """
    Predict oriC using nucleotide skew analysis + BLAST against DoriC database
    
    Args:
        sequence: Input DNA sequence
        annotation_df: Existing annotations dataframe
        doric_fasta: Path to DoriC database FASTA
        min_length: Minimum length for BLAST hits (default: 100bp)
        location_buffer: Buffer between multiple hits (default: 50bp)
        window_size: Window size for skew analysis (default: 500bp)
        sequence_id: Sequence identifier
        max_search_region: Maximum region around predicted center for BLAST (default: 3000bp)
    
    Returns:
        Updated annotation dataframe with oriC predictions
    """
    print("🔍 Predicting oriC using skew analysis + DoriC database...")
    
    seq = sequence.lower()
    
    # === STEP 1: Nucleotide Skew Analysis ===
    disparity_curves = {
        "GC": {"positive": "g", "negative": "c"},
        "AT": {"positive": "a", "negative": "t"},
        "MK": {"positive": "a", "negative": "c"},
        "RY": {"positive": "a", "negative": "g"},
    }

    results = {}
    for name, groups in disparity_curves.items():
        curve = []
        pos = groups["positive"]
        neg = groups["negative"]
        for i in range(0, len(seq) - window_size + 1, window_size):
            window = seq[i:i + window_size]
            disparity = window.count(pos) - window.count(neg)
            curve.append(disparity)
        results[name] = np.cumsum(curve)

    # Normalize curves and find minimum
    normalized_curves = {k: (v - np.mean(v)) / np.std(v) for k, v in results.items()}
    combined_score = sum(normalized_curves.values())
    min_index = np.argmin(combined_score)
    predicted_center = (min_index * window_size + (min_index + 1) * window_size) // 2
    
    print(f"   📊 Skew analysis predicts oriC center at position: {predicted_center}")

    # === STEP 2: Extract Search Region for BLAST ===
    half_search = max_search_region // 2
    search_start = max(0, predicted_center - half_search)
    search_end = min(len(seq), predicted_center + half_search)
    search_seq = seq[search_start:search_end]
    
    print(f"   🔍 Searching {len(search_seq)}bp region ({search_start}-{search_end}) against DoriC database")

    # === STEP 3: Create temporary FASTA file for BLAST ===
    temp_fasta = f"{sequence_id}_oriC_search.fasta"
    with open(temp_fasta, "w") as f:
        f.write(f">{sequence_id}_search\n{search_seq}\n")

    # === STEP 4: Set up BLAST database ===
    db_folder = Path("oric_db_folder")
    db_folder.mkdir(exist_ok=True)
    doric_db_name = str(db_folder / "doric_temp_db")

    # === STEP 5: Run BLAST ===
    blast_output = "oric_blast_hits.csv"
    try:
        os.system(f"makeblastdb -in {doric_fasta} -dbtype nucl -out {doric_db_name}")
        os.system(
            f"blastn -query {temp_fasta} -db {doric_db_name} "
            f"-evalue 1e-3 -outfmt '10 qseqid sseqid pident length mismatch gapopen "
            f"qstart qend sstart send evalue bitscore stitle' "
            f"-max_target_seqs 10 -out {blast_output}"
        )
    except Exception as e:
        print(f"❌ BLAST error: {e}")

    # === STEP 6: Parse BLAST results ===
    header = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    
    oriC_predictions = []
    
    if os.path.exists(blast_output):
        try:
            doric_df = pd.read_csv(blast_output, header=None, names=header)
            
            if not doric_df.empty:
                # Filter by minimum length and sort by bit score
                doric_df = doric_df[doric_df["length"] >= min_length].sort_values("bitscore", ascending=False)
                
                print(f"   📊 Found {len(doric_df)} DoriC BLAST hits (≥{min_length}bp)")
                
                # Process each significant hit
                filtered_hits = []
                last_end = -float("inf")
                
                for _, row in doric_df.iterrows():
                    # Convert coordinates back to original sequence
                    hit_start = search_start + min(row.qstart, row.qend)
                    hit_end = search_start + max(row.qstart, row.qend)
                    
                    # Check if this hit is far enough from previous hits
                    if hit_start - last_end > location_buffer:
                        filtered_hits.append({
                            'start': hit_start,
                            'end': hit_end,
                            'strand': '+' if row.qstart <= row.qend else '-',
                            'identity': row.pident,
                            'evalue': row.evalue,
                            'hit_description': row.stitle,
                            'length': row.length
                        })
                        last_end = hit_end
                        
                        print(f"   ✅ DoriC hit: {hit_start}-{hit_end} ({row.length}bp, {row.pident:.1f}% identity)")
                
                oriC_predictions = filtered_hits
            else:
                print("   ⚠️ No significant DoriC BLAST hits found")
                
        except Exception as e:
            print(f"   ❌ Error parsing BLAST results: {e}")

    # === STEP 7: Create oriC annotations ===
    if oriC_predictions:
        # Use BLAST hits
        for i, hit in enumerate(oriC_predictions):
            oriC_seq = sequence[hit['start']-1:hit['end']]  # Convert to 0-based indexing
            if hit['strand'] == '-':
                oriC_seq = str(Seq(oriC_seq).reverse_complement())
            
            oriC_entry = {
                "Gene Name": f"oriV_{i+1}" if len(oriC_predictions) > 1 else "oriV",
                "Product": f"Origin of replication (DoriC: {hit['identity']:.1f}% identity)",
                "Start": hit['start'],
                "End": hit['end'],
                "Strand": hit['strand'],
                "Category": "Origin of Replication",
                "Translation": oriC_seq.upper(),
                "feature type": "misc_feature"
            }
            
            # Add missing columns
            for col in annotation_df.columns:
                if col not in oriC_entry:
                    oriC_entry[col] = ""
            
            annotation_df = pd.concat([annotation_df, pd.DataFrame([oriC_entry])], ignore_index=True)
            
    else:
        # Fallback: Use predicted center with small region
        print("   🔄 No DoriC hits found, using skew-predicted region")
        
        # Check for DnaA motifs in smaller region around predicted center
        fallback_size = 300  # Much smaller fallback region
        fallback_start = max(0, predicted_center - fallback_size//2)
        fallback_end = min(len(seq), predicted_center + fallback_size//2)
        fallback_seq = seq[fallback_start:fallback_end]
        
        # Look for DnaA motifs
        dnaa_motifs = [
            "ttatccaca", "ttatccaaa", "ttgtccaca", "ttatgcaca", 
            "ttttccaca", "ttatccata", "ttatccacc"
        ]
        motif_found = any(re.search(m, fallback_seq) for m in dnaa_motifs)
        
        confidence = "high" if motif_found else "low"
        product_desc = f"Predicted origin of replication (skew analysis, {confidence} confidence)"
        
        oriC_entry = {
            "Gene Name": "oriV",
            "Product": product_desc,
            "Start": fallback_start + 1,  # Convert to 1-based
            "End": fallback_end,
            "Strand": "+",
            "Category": "Origin of Replication",
            "Translation": fallback_seq.upper(),
            "feature type": "misc_feature"
        }
        
        # Add missing columns
        for col in annotation_df.columns:
            if col not in oriC_entry:
                oriC_entry[col] = ""

        annotation_df = pd.concat([annotation_df, pd.DataFrame([oriC_entry])], ignore_index=True)
        
        print(f"   📍 Predicted oriV: {fallback_start + 1}-{fallback_end} ({fallback_size}bp, DnaA motifs: {motif_found})")

    # === STEP 8: Clean up ===
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
    if os.path.exists(blast_output):
        os.remove(blast_output)
    if db_folder.exists():
        shutil.rmtree(db_folder)

    num_predictions = len(oriC_predictions) if oriC_predictions else 1
    print(f"✅ oriC prediction complete. Added {num_predictions} origin(s) of replication")
    
    return annotation_df



def predict_oriT(sequence,annotation_df,oriT_fastas,min_identity=80.0,min_length=100,location_buffer=50):
    os.makedirs("orit_db_folder", exist_ok=True)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w") as tmp_fasta:
        tmp_fasta.write(">query_seq\n")
        tmp_fasta.write(sequence)
        tmp_fasta_path = tmp_fasta.name

    all_hits = []
    blast_output_files = []

    for oriT_fasta in oriT_fastas:
        db_name = Path(oriT_fasta).stem
        db_prefix = os.path.join("orit_db_folder", db_name)
        out_file = db_prefix + "_hits.csv"
        blast_output_files.append(out_file)

        os.system(f"makeblastdb -in {oriT_fasta} -dbtype nucl -out {db_prefix}")
        os.system(
            f"blastn -query {tmp_fasta_path} -db {db_prefix} "
            f"-evalue 1e-5 -outfmt '10 qseqid sseqid pident length mismatch gapopen "
            f"qstart qend sstart send evalue bitscore stitle' "
            f"-max_target_seqs 10 -out {out_file}"
        )

        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
        ]
        try:
            df = pd.read_csv(out_file, header=None, names=columns)
            df = df[(df["pident"] >= min_identity) & (df["length"] >= min_length)]
            df["source_db"] = Path(oriT_fasta).name
            all_hits.append(df)
        except Exception:
            continue

    # Cleanup temporary query and output files
    os.remove(tmp_fasta_path)
    for f in blast_output_files:
        if os.path.exists(f):
            os.remove(f)

    # Cleanup the entire orit_db_folder
    shutil.rmtree("orit_db_folder", ignore_errors=True)

    if not all_hits:
        print("[INFO] No oriT hits found.")
        return annotation_df

    all_hits_df = pd.concat(all_hits).sort_values("qstart")

    filtered = []
    last_end = -float("inf")
    for _, row in all_hits_df.iterrows():
        if row.qstart - last_end > location_buffer:
            filtered.append(row)
            last_end = row.qend
    hits_df = pd.DataFrame(filtered)

    for _, hit in hits_df.iterrows():
        start = min(hit.qstart, hit.qend) - 1
        end = max(hit.qstart, hit.qend)
        seq = sequence[start:end]
        if hit.qstart > hit.qend:
            seq = str(Seq(seq).reverse_complement())

        new_row = {
            "Gene Name": "oriT",
            "Product": f"Predicted origin of transfer ({hit['source_db']})",
            "Start": start + 1,
            "End": end,
            "Strand": "+" if hit.qstart <= hit.qend else "-",
            "Category": "Transfer",
            "Translation": seq
        }
        for col in annotation_df.columns:
            if col not in new_row:
                new_row[col] = ""
        annotation_df = pd.concat([annotation_df, pd.DataFrame([new_row])], ignore_index=True)

    print(f"[INFO] oriT prediction complete. Appended {len(hits_df)} oriT regions.")
    return annotation_df


def predict_replicons(
    sequence,
    annotation_df,
    replicon_fasta,
    min_identity=60.0,
    min_length=50,
    location_buffer=50
):
    # Create temp query FASTA
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w") as tmp_fasta:
        tmp_fasta.write(">query_seq\n")
        tmp_fasta.write(sequence)
        tmp_fasta_path = tmp_fasta.name

    # Create BLAST database folder
    os.makedirs("replicon_db_folder", exist_ok=True)
    db_prefix = os.path.join("replicon_db_folder", Path(replicon_fasta).stem)
    blast_output = db_prefix + "_hits.csv"

    # Create BLAST DB
    os.system(f"makeblastdb -in {replicon_fasta} -dbtype nucl -out {db_prefix}")

    # Run BLAST
    blast_cmd = (
        f"blastn -query {tmp_fasta_path} -db {db_prefix} "
        f"-evalue 1e-5 -outfmt '10 qseqid sseqid pident length mismatch gapopen "
        f"qstart qend sstart send evalue bitscore stitle' "
        f"-max_target_seqs 10 -out {blast_output}"
    )
    os.system(blast_cmd)

    # Load and filter BLAST hits
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    try:
        df = pd.read_csv(blast_output, header=None, names=columns)
        df = df[(df["pident"] >= min_identity) & (df["length"] >= min_length)]
    except Exception as e:
        print(f"[ERROR] Could not read BLAST hits: {e}")
        df = pd.DataFrame()

    # Clean up temp files
    os.remove(tmp_fasta_path)
    if os.path.exists(blast_output):
        os.remove(blast_output)
    shutil.rmtree("replicon_db_folder", ignore_errors=True)

    if df.empty:
        print("[INFO] No replicon hits found.")
        return annotation_df

    # Filter overlapping hits
    df = df.sort_values("qstart")
    filtered = []
    last_end = -float("inf")
    for _, row in df.iterrows():
        if row.qstart - last_end > location_buffer:
            filtered.append(row)
            last_end = row.qend
    hits_df = pd.DataFrame(filtered)

    # Add new replicon entries
    for _, hit in hits_df.iterrows():
        start = min(hit.qstart, hit.qend) - 1
        end = max(hit.qstart, hit.qend)
        seq = sequence[start:end]
        if hit.qstart > hit.qend:
            seq = str(Seq(seq).reverse_complement())

        new_row = {
            "Gene Name": f"Replicon_{hit['sseqid']}",
            "Product": f"Predicted plasmid replicon ({hit['sseqid']})",
            "Start": start + 1,
            "End": end,
            "Strand": "+" if hit.qstart <= hit.qend else "-",
            "Category": "Replication",
            "Translation": seq
        }
        for col in annotation_df.columns:
            if col not in new_row:
                new_row[col] = ""
        annotation_df = pd.concat([annotation_df, pd.DataFrame([new_row])], ignore_index=True)

    print(f"[INFO] Replicon prediction complete. Appended {len(hits_df)} replicon regions.")
    return annotation_df


def predict_transposons(sequence,annotation_df,transposon_fastas,min_identity=80.0,min_length=100,location_buffer=50):
    # Write query sequence to temporary FASTA
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w") as tmp_fasta:
        tmp_fasta.write(">query_seq\n")
        tmp_fasta.write(sequence)
        tmp_fasta_path = tmp_fasta.name

    # Setup DB folder
    db_folder = "transposon_db_folder"
    os.makedirs(db_folder, exist_ok=True)
    blast_output_files = []
    all_hits = []

    for transposon_fasta in transposon_fastas:
        db_prefix = os.path.join(db_folder, Path(transposon_fasta).stem)
        out_file = db_prefix + "_hits.csv"
        blast_output_files.append(out_file)

        # Make BLAST DB and run BLAST
        os.system(f"makeblastdb -in {transposon_fasta} -dbtype nucl -out {db_prefix}")
        os.system(
            f"blastn -query {tmp_fasta_path} -db {db_prefix} "
            f"-evalue 1e-5 -outfmt '10 qseqid sseqid pident length mismatch gapopen "
            f"qstart qend sstart send evalue bitscore stitle' "
            f"-max_target_seqs 10 -out {out_file}"
        )

        # Load BLAST results
        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
        ]
        try:
            df = pd.read_csv(out_file, header=None, names=columns)
            df = df[(df["pident"] >= min_identity) & (df["length"] >= min_length)]
            df["source_db"] = Path(transposon_fasta).name
            all_hits.append(df)
        except Exception as e:
            print(f"[WARNING] Could not process {transposon_fasta}: {e}")

    # Clean up temporary query and results
    os.remove(tmp_fasta_path)
    for f in blast_output_files:
        if os.path.exists(f):
            os.remove(f)
    shutil.rmtree(db_folder, ignore_errors=True)

    if not all_hits:
        print("[INFO] No transposon hits found.")
        return annotation_df

    all_hits_df = pd.concat(all_hits).sort_values("qstart")

    # Remove overlapping hits
    filtered = []
    last_end = -float("inf")
    for _, row in all_hits_df.iterrows():
        if row.qstart - last_end > location_buffer:
            filtered.append(row)
            last_end = row.qend
    hits_df = pd.DataFrame(filtered)

    # Add hits to annotation DataFrame
    for _, hit in hits_df.iterrows():
        start = min(hit.qstart, hit.qend) - 1
        end = max(hit.qstart, hit.qend)
        seq = sequence[start:end]
        if hit.qstart > hit.qend:
            seq = str(Seq(seq).reverse_complement())

        name_part = hit['sseqid'].split("-")[0].split("_")[-1]
        is_element = name_part.startswith("IS")

        new_row = {
            "Gene Name": name_part,
            "Product": "Predicted insertion sequence" if is_element else "Predicted transposon",
            "Start": start + 1,
            "End": end,
            "Strand": "+" if hit.qstart <= hit.qend else "-",
            "Category": "Mobile Element",
            "Translation": seq
        }
        for col in annotation_df.columns:
            if col not in new_row:
                new_row[col] = ""
        annotation_df = pd.concat([annotation_df, pd.DataFrame([new_row])], ignore_index=True)

    print(f"[INFO] Transposon prediction complete. Appended {len(hits_df)} transposon regions.")
    return annotation_df


def run_infernal_on_sequence(sequence, rfam_cm_path):
    # Check for Infernal tools
    if not shutil.which("cmscan") or not shutil.which("cmpress"):
        raise EnvironmentError("❌ Infernal (cmscan/cmpress) is not installed or not in PATH.")

    # Auto-run cmpress if index files are missing
    if not all(os.path.exists(f"{rfam_cm_path}.{ext}") for ext in ["i1f", "i1i", "i1m", "i1p"]):
        print("📦 Running cmpress to index Rfam.cm...")
        subprocess.run(["cmpress", rfam_cm_path], check=True)

    # Use a temp directory to store the sequence FASTA and cmscan output
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "query.fasta")
        cmscan_out = os.path.join(tmpdir, "cmscan.out")

        with open(fasta_path, "w") as f:
            f.write(">query\n")
            f.write(sequence + "\n")

        with open(cmscan_out, "w") as output_file:
            subprocess.run([
                "cmscan", "--cut_ga", "--rfam", "--nohmmonly",
                rfam_cm_path, fasta_path
            ], stdout=output_file, check=True)


        # Parse the output
        with open(cmscan_out, "r") as f:
            lines = f.readlines()

        records = []
        for i, line in enumerate(lines):
            if line.startswith(">>"):
                parts = line.strip().split(None, 2)
                if len(parts) == 3:
                    family, description = parts[1], parts[2]
                    for subline in lines[i+1:]:
                        if subline.startswith("  ("):
                            match = re.search(r'(\d+)\s+(\d+)\s+([+-])', subline)
                            if match:
                                start, end, strand = match.groups()
                                records.append([family, description, int(start), int(end), strand])
                            break

        return pd.DataFrame(records, columns=["ncRNA Family", "Description", "Start", "End", "Strand"])
    

import pandas as pd
import math
import re
from Bio import SeqIO
from Bio.Seq import Seq

# --- PWM DEFINITIONS ---
PWM_35 = {
    0: {'T': 1.0}, 1: {'T': 1.0}, 2: {'G': 0.7, 'A': 0.3},
    3: {'A': 0.6, 'T': 0.4}, 4: {'C': 1.0}, 5: {'A': 1.0}
}

PWM_10 = {
    0: {'T': 1.0}, 1: {'A': 1.0}, 2: {'T': 1.0},
    3: {'A': 1.0}, 4: {'A': 0.7, 'T': 0.3}, 5: {'T': 1.0}
}

# --- UTILITIES ---
def score_pwm(seq, pwm):
    score = 0
    for i, base in enumerate(seq.upper()):
        probs = pwm.get(i, {})
        p = probs.get(base, 0.01)
        score += math.log2(p / 0.25)
    return round(score, 2)

def hamming(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# --- PROMOTER DETECTION ---
def detect_promoters(sequence, gene_df, upstream_bp=100, max_mismatches=2):
    sequence = sequence.upper()
    promoters = []

    for _, row in gene_df.iterrows():
        gene_name = row["Gene Name"]
        product = row["Product"]
        category = row["Category"]
        strand = row["Strand"]

        if strand == "+":
            start = max(0, row["Start"] - upstream_bp)
            region = sequence[start:row["Start"]]
            #print(f"[{gene_name}] Strand {strand} — Region Length: {len(region)} — From {start} to {row['Start']}")
            abs_offset = start
        elif strand == "-":
            end = min(len(sequence), row["End"] + upstream_bp)
            region = str(Seq(sequence[row["End"]:end]).reverse_complement())
            #print(f"[{gene_name}] Strand {strand} — Region Length: {len(region)} — From {start} to {row['Start']}")
            abs_offset = row["End"]
        else:
            continue

        for i in range(len(region) - 50):
            box_35 = region[i:i+6]
            if hamming(box_35, "TTGACA") <= max_mismatches:
                for spacer in range(15, 20):
                    j = i + 6 + spacer
                    box_10 = region[j:j+6]
                    if len(box_10) < 6:
                        continue
                    if hamming(box_10, "TATAAT") <= max_mismatches:
                        score35 = score_pwm(box_35, PWM_35)
                        score10 = score_pwm(box_10, PWM_10)
                        total_score = score35 + score10
                        promoters.append({
                            "Gene": gene_name,
                            "Product": product,
                            "Category": category,
                            "Gene Strand": strand,
                            "Promoter Region Start": abs_offset + i,
                            "Promoter Region End": abs_offset + j + 6,
                            "-35 Box": box_35,
                            "-10 Box": box_10,
                            "Spacer": spacer,
                            "Score -35": score35,
                            "Score -10": score10,
                            "Total Score": total_score
                        })
    promoter_df = pd.DataFrame(promoters)
    filtered = promoter_df[
        (promoter_df["Total Score"] >= 4) &
        (promoter_df["Spacer"].between(15, 19)) &
        (promoter_df["Score -35"] >= 1.5) &
        (promoter_df["Score -10"] >= 1.5)
    ]
    return filtered


# ========================================================================================
# UNIPROT BLAST ANNOTATION FUNCTIONS
# ========================================================================================

def prepare_uniprot_blast_database(uniprot_tsv_path, output_dir="uniprot_blast_db"):
    """
    Prepare BLAST database from UniProt TSV file
    
    Args:
        uniprot_tsv_path (str): Path to UniProt TSV file
        output_dir (str): Directory to store BLAST database
    
    Returns:
        str: Path to BLAST database prefix
    """
    print("🔄 Preparing UniProt BLAST database...")
    
    # Create output directory
    db_dir = Path(output_dir)
    db_dir.mkdir(exist_ok=True)
    
    fasta_path = db_dir / "uniprot_sequences.fasta"
    db_prefix = db_dir / "uniprot_db"
    
    # Check if database already exists
    if (db_dir / "uniprot_db.pin").exists():
        print("✅ UniProt BLAST database already exists")
        return str(db_prefix)
    
    # Read UniProt TSV and create FASTA
    print("📖 Reading UniProt TSV file...")
    df = pd.read_csv(uniprot_tsv_path, sep='\t')
    
    print(f"📊 Found {len(df)} UniProt entries")
    
    # Create FASTA file from UniProt sequences
    with open(fasta_path, 'w') as fasta_out:
        for idx, row in df.iterrows():
            entry_id = row.get('Entry', f'UNK_{idx}')
            sequence = row.get('Sequence', '')
            
            if pd.notnull(sequence) and sequence.strip():
                # Header includes entry ID for later lookup
                fasta_out.write(f">{entry_id}\n{sequence.strip()}\n")
    
    print(f"✅ Created FASTA with {len(df)} sequences")
    
    # Build BLAST database
    print("🛠️ Building UniProt BLAST database...")
    try:
        subprocess.run([
            "makeblastdb",
            "-in", str(fasta_path),
            "-dbtype", "prot",
            "-out", str(db_prefix),
            "-title", "UniProt_Plasmid_Database"
        ], check=True, capture_output=True, text=True)
        
        print("✅ UniProt BLAST database created successfully")
        return str(db_prefix)
        
    except subprocess.CalledProcessError as e:
        raise Exception(f"Failed to create UniProt BLAST database: {e.stderr}")
    except FileNotFoundError:
        raise Exception("makeblastdb not found! Please install BLAST+ tools.")

def run_blastp_against_uniprot(query_fasta, uniprot_db_prefix, output_file, max_target_seqs=5, evalue=1e-5):
    """
    Run BLASTP search against UniProt database
    
    Args:
        query_fasta (str): Path to query FASTA file
        uniprot_db_prefix (str): UniProt BLAST database prefix
        output_file (str): Output file for BLAST results
        max_target_seqs (int): Maximum number of target sequences
        evalue (float): E-value threshold
    """
    print(f"🔍 Running BLASTP against UniProt database...")
    
    try:
        subprocess.run([
            "blastp",
            "-query", query_fasta,
            "-db", uniprot_db_prefix,
            "-out", output_file,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-max_target_seqs", str(max_target_seqs),
            "-evalue", str(evalue),
            "-num_threads", "4"
        ], check=True, capture_output=True, text=True)
        
        print("✅ BLASTP search completed")
        
    except subprocess.CalledProcessError as e:
        raise Exception(f"BLASTP failed: {e.stderr}")
    except FileNotFoundError:
        raise Exception("blastp not found! Please install BLAST+ tools.")

def parse_uniprot_blast_results(blast_output_file, uniprot_tsv_path, min_identity=30, min_coverage=50):
    """
    Parse BLAST results and create mapping to UniProt annotations
    
    Args:
        blast_output_file (str): BLAST output file
        uniprot_tsv_path (str): UniProt TSV file for annotation lookup
        min_identity (float): Minimum identity percentage
        min_coverage (float): Minimum query coverage percentage
    
    Returns:
        dict: Mapping of query IDs to UniProt annotations
    """
    print("📊 Parsing BLAST results...")
    
    # Read UniProt annotations
    uniprot_df = pd.read_csv(uniprot_tsv_path, sep='\t')
    uniprot_dict = {}
    
    # Create lookup dictionary
    for idx, row in uniprot_df.iterrows():
        entry_id = row.get('Entry', '')
        if entry_id:
            uniprot_dict[entry_id] = {
                'protein_names': row.get('Protein names', ''),
                'gene_names': row.get('Gene Names', ''),
                'organism': row.get('Organism', ''),
                'function': row.get('Function [CC]', ''),
                'go_biological': row.get('Gene Ontology (biological process)', ''),
                'go_cellular': row.get('Gene Ontology (cellular component)', ''),
                'go_molecular': row.get('Gene Ontology (molecular function)', ''),
                'pubmed_id': row.get('PubMed ID', '')
            }
    
    # Parse BLAST results
    blast_annotations = {}
    
    if not os.path.exists(blast_output_file):
        print("⚠️ No BLAST results file found")
        return blast_annotations
    
    with open(blast_output_file, 'r') as blast_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query_id = fields[0]
                subject_id = fields[1]
                identity = float(fields[2])
                length = int(fields[3])
                qstart = int(fields[6])
                qend = int(fields[7])
                evalue = float(fields[10])
                bitscore = float(fields[11])
                
                # Apply filters
                if identity >= min_identity:
                    # Only keep best hit per query
                    if query_id not in blast_annotations:
                        if subject_id in uniprot_dict:
                            annotation = uniprot_dict[subject_id].copy()
                            annotation.update({
                                'blast_identity': identity,
                                'blast_evalue': evalue,
                                'blast_bitscore': bitscore,
                                'uniprot_entry': subject_id
                            })
                            blast_annotations[query_id] = annotation
    
    print(f"✅ Found annotations for {len(blast_annotations)} ORFs")
    return blast_annotations

def extract_protein_sequences_from_dataframe(dataframe):
    """
    Extract protein sequences from the dataframe and create FASTA file
    
    Args:
        dataframe (pd.DataFrame): DataFrame with ORF information
    
    Returns:
        str: Path to temporary FASTA file
    """
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    
    protein_count = 0
    for idx, row in dataframe.iterrows():
        if row.get('Type', '') == 'CDS':
            sequence = row.get('Sequence', '')
            if sequence:
                # Translate DNA to protein if needed
                try:
                    # Check if it's already protein or needs translation
                    if re.match(r'^[ACGT]+$', sequence.upper()):
                        # It's DNA, translate it
                        protein_seq = str(Seq(sequence).translate().rstrip('*'))
                    else:
                        # Assume it's already protein
                        protein_seq = sequence.rstrip('*')
                    
                    if len(protein_seq) > 10:  # Minimum protein length
                        orf_id = f"ORF_{idx}_{row.get('Start', 'UNK')}_{row.get('End', 'UNK')}"
                        temp_fasta.write(f">{orf_id}\n{protein_seq}\n")
                        protein_count += 1
                        
                except Exception as e:
                    print(f"⚠️ Could not process sequence for ORF {idx}: {e}")
                    continue
    
    temp_fasta.close()
    print(f"📄 Created query FASTA with {protein_count} protein sequences")
    return temp_fasta.name

def uniprot_blast(final_dataframe_fixed, uniprot_tsv_path, min_identity=30, evalue=1e-5):
    """
    Main function to perform UniProt BLAST annotation
    
    Args:
        final_dataframe_fixed (pd.DataFrame): DataFrame with ORF annotations
        uniprot_tsv_path (str): Path to UniProt TSV file
        min_identity (float): Minimum identity percentage for accepting hits
        evalue (float): E-value threshold for BLAST
    
    Returns:
        pd.DataFrame: Updated dataframe with UniProt annotations
    """
    print("🧬 Starting UniProt BLAST annotation...")
    print("⚠️ This may take several minutes depending on the number of ORFs...")
    
    # Create a copy of the dataframe
    updated_df = final_dataframe_fixed.copy()
    
    # Step 1: Prepare UniProt BLAST database
    try:
        uniprot_db_prefix = prepare_uniprot_blast_database(uniprot_tsv_path)
    except Exception as e:
        print(f"❌ Failed to prepare UniProt database: {e}")
        return updated_df
    
    # Step 2: Extract protein sequences from ORFs
    try:
        query_fasta_path = extract_protein_sequences_from_dataframe(updated_df)
    except Exception as e:
        print(f"❌ Failed to extract protein sequences: {e}")
        return updated_df
    
    # Step 3: Run BLAST search
    blast_output_file = "uniprot_blast_results.txt"
    try:
        run_blastp_against_uniprot(query_fasta_path, uniprot_db_prefix, blast_output_file, evalue=evalue)
    except Exception as e:
        print(f"❌ BLAST search failed: {e}")
        # Clean up
        if os.path.exists(query_fasta_path):
            os.unlink(query_fasta_path)
        return updated_df
    
    # Step 4: Parse BLAST results and get annotations
    try:
        blast_annotations = parse_uniprot_blast_results(blast_output_file, uniprot_tsv_path, min_identity)
    except Exception as e:
        print(f"❌ Failed to parse BLAST results: {e}")
        blast_annotations = {}
    
    # Step 5: Update dataframe with UniProt annotations
    updated_count = 0
    for idx, row in updated_df.iterrows():
        if row.get('Type', '') == 'CDS':
            orf_id = f"ORF_{idx}_{row.get('Start', 'UNK')}_{row.get('End', 'UNK')}"
            
            if orf_id in blast_annotations:
                annotation = blast_annotations[orf_id]
                
                # Update Gene Name (use first gene name if multiple)
                gene_names = annotation.get('gene_names', '')
                if gene_names and pd.notnull(gene_names):
                    # Take first gene name if multiple (space or comma separated)
                    first_gene = re.split(r'[,\s]+', str(gene_names))[0].strip()
                    if first_gene:
                        updated_df.at[idx, 'Gene Name'] = first_gene
                
                # Update Product (use protein names)
                protein_names = annotation.get('protein_names', '')
                if protein_names and pd.notnull(protein_names):
                    # Clean up protein names (remove EC numbers, etc.)
                    clean_protein_name = re.sub(r'\s*\(EC[^)]*\)', '', str(protein_names))
                    clean_protein_name = clean_protein_name.split(';')[0].strip()  # Take first name
                    if clean_protein_name:
                        updated_df.at[idx, 'Product'] = clean_protein_name
                
                # Add additional UniProt information as new columns (optional)
                updated_df.at[idx, 'UniProt_Entry'] = annotation.get('uniprot_entry', '')
                updated_df.at[idx, 'UniProt_Organism'] = annotation.get('organism', '')
                updated_df.at[idx, 'BLAST_Identity'] = annotation.get('blast_identity', '')
                updated_df.at[idx, 'BLAST_Evalue'] = annotation.get('blast_evalue', '')
                
                updated_count += 1
    
    print(f"✅ Updated {updated_count} ORFs with UniProt annotations")
    
    # Clean up temporary files
    if os.path.exists(query_fasta_path):
        os.unlink(query_fasta_path)
    if os.path.exists(blast_output_file):
        os.unlink(blast_output_file)
    
    print("🎉 UniProt BLAST annotation completed!")
    return updated_df

def add_uniprot_annotation_option(final_dataframe_fixed, uniprot_tsv_path=None, 
                                  run_uniprot_blast=False, min_identity=30):
    """
    Wrapper function to optionally add UniProt annotations
    
    Args:
        final_dataframe_fixed (pd.DataFrame): DataFrame with ORF annotations
        uniprot_tsv_path (str): Path to UniProt TSV file
        run_uniprot_blast (bool): Whether to run UniProt BLAST
        min_identity (float): Minimum identity percentage
    
    Returns:
        pd.DataFrame: Potentially updated dataframe
    """
    if run_uniprot_blast and uniprot_tsv_path:
        if not os.path.exists(uniprot_tsv_path):
            print(f"❌ UniProt TSV file not found: {uniprot_tsv_path}")
            return final_dataframe_fixed
        
        print("🚀 Running optional UniProt BLAST annotation...")
        return uniprot_blast(final_dataframe_fixed, uniprot_tsv_path, min_identity)
    else:
        print("⏭️ Skipping UniProt BLAST annotation")
        return final_dataframe_fixed
    
    # ========================================================================================
# ========================================================================================
# REPLACE ALL THE INTERGENIC FUNCTIONS IN blast_search.py WITH THESE SIMPLE VERSIONS
# ========================================================================================
# ========================================================================================
# REPLACE ALL INTERGENIC FUNCTIONS WITH THESE CLEAN SINGLE-BLAST VERSIONS
# ========================================================================================

def find_intergenic_regions_simple(annotation_df, sequence, min_gap_size=100):
    """Find intergenic regions using simple logic"""
    # Get CDS features and sort by start position
    cds_features = annotation_df[annotation_df['feature type'] == 'CDS'].copy()
    cds_features = cds_features.sort_values('Start').reset_index(drop=True)
    
    if len(cds_features) < 2:
        return []
    
    intergenic_regions = []
    
    for i in range(len(cds_features) - 1):
        prev_cds = cds_features.iloc[i]
        next_cds = cds_features.iloc[i + 1]
        
        gap_start = prev_cds['End']
        gap_end = next_cds['Start']
        gap_size = gap_end - gap_start
        
        if gap_size >= min_gap_size:
            # Extract sequence
            gap_sequence = sequence[gap_start-1:gap_end-1]  # Convert to 0-based indexing
            
            intergenic_regions.append({
                'region_id': f"intergenic_{i+1}",
                'start': gap_start,
                'end': gap_end,
                'length': gap_size,
                'sequence': gap_sequence
            })
    
    return intergenic_regions


def blast_intergenic_regions_single(intergenic_regions, blast_db_prefix, temp_dir, max_workers=4):
    """BLAST intergenic regions - SINGLE BLAST per region"""
    print("🔍 Running BLAST on intergenic regions (single BLAST per region)...")
    
    if not intergenic_regions:
        return []
    
    tmp_dir = Path(temp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    
    blast_header = ",".join([
        "query_id", "subject_id", "identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end",
        "s_start", "s_end", "evalue", "bit_score"
    ])
    
    # Create query files and run BLAST - ONE per region
    tasks = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for region in intergenic_regions:
            # Submit BLAST job for region (BLASTX handles both strands automatically)
            tasks.append(executor.submit(
                run_blast_for_single_intergenic_region, 
                region['region_id'], region['sequence'], 
                str(tmp_dir), blast_db_prefix, blast_header
            ))
        
        # Wait for all tasks to complete
        for f in as_completed(tasks):
            try:
                f.result()
            except Exception as e:
                print(f"❌ Error in BLAST task: {e}")
    
    print("✅ BLAST searches completed")
    return tmp_dir


def run_blast_for_single_intergenic_region(region_id, sequence, tmp_dir, blast_db_prefix, blast_header):
    """Run BLAST for a single intergenic region - BLASTX handles both strands"""
    tmp_dir = Path(tmp_dir)
    query_file = tmp_dir / f"{region_id}.fasta"
    result_file = tmp_dir / f"{region_id}_blast.csv"
    
    # Write query FASTA (only forward sequence - BLASTX does 6-frame translation)
    with open(query_file, "w") as f:
        f.write(f">{region_id}\n{sequence}\n")
    
    # Run BLASTX (automatically searches both strands)
    blast_cmd = (
        f"blastx -query {query_file} -db {blast_db_prefix} "
        f"-out {result_file} -outfmt '10' -max_target_seqs 5 -evalue 1e-5"
    )
    os.system(blast_cmd)
    
    # Add header
    if result_file.exists() and result_file.stat().st_size > 0:
        with open(result_file, "r+") as f:
            content = f.read()
            f.seek(0)
            f.write(blast_header + "\n" + content)


def parse_intergenic_blast_single(blast_results_dir, database_csv_path, intergenic_regions, sequence):
    """Parse BLAST results - SINGLE BLAST version - SILENT"""
    blast_results_dir = Path(blast_results_dir)
    
    # Load database for annotations
    try:
        db_df = pd.read_csv(database_csv_path)
    except Exception as e:
        return pd.DataFrame()
    
    all_hits = []
    
    # Process each intergenic region
    for region in intergenic_regions:
        region_id = region['region_id']
        blast_file = blast_results_dir / f"{region_id}_blast.csv"
        
        if not blast_file.exists():
            continue
        
        try:
            blast_df = pd.read_csv(blast_file)
            if blast_df.empty:
                continue
                
            # Filter by identity
            blast_df = blast_df[blast_df['identity'] >= 40]
            if blast_df.empty:
                continue
            
            # Get best hit (highest bit score)
            best_hit = blast_df.loc[blast_df['bit_score'].idxmax()]
            
            # Get gene info from database
            try:
                gene_info = db_df.iloc[int(best_hit['subject_id'])]
                gene_name = gene_info.get('Gene Name', f"intergenic_{region_id}")
                product = gene_info.get('Product', 'Predicted protein')
                category = gene_info.get('Category', 'Unknown')
            except:
                gene_name = f"intergenic_{region_id}"
                product = 'Predicted protein'
                category = 'Unknown'
            
            # Simple coordinate calculation
            hit_start = best_hit['q_start']
            hit_end = best_hit['q_end']
            
            # Calculate genomic coordinates
            genomic_start = region['start'] + min(hit_start, hit_end) - 1
            genomic_end = region['start'] + max(hit_start, hit_end) - 1
            
            # Proper strand determination from BLAST coordinates
            if hit_start <= hit_end:
                strand = '+'  # Forward orientation
            else:
                strand = '-'  # Reverse orientation
            
            # Extract sequence
            hit_sequence = sequence[genomic_start-1:genomic_end]
            if strand == '-':
                hit_sequence = str(Seq(hit_sequence).reverse_complement())
            
            all_hits.append({
                'Gene Name': gene_name,
                'Product': product,
                'Start': genomic_start,
                'End': genomic_end,
                'Strand': strand,
                'Category': category,
                'Translation': hit_sequence,
                'feature type': 'CDS',
                'Intergenic_Flag': 1
            })
            
        except Exception as e:
            # Silent error handling - no more error messages
            continue
    
    return pd.DataFrame(all_hits)

def filter_overlapping_intergenic_genes(combined_df, min_length=90, overlap_threshold=50):
    """Filter out intergenic genes that overlap with original CDS genes or are too short - SILENT"""
    # Separate original and intergenic genes
    original_genes = combined_df[combined_df['Intergenic_Flag'] == 0].copy()
    intergenic_genes = combined_df[combined_df['Intergenic_Flag'] == 1].copy()
    
    # Filter 1: Remove short intergenic genes
    intergenic_genes['Length'] = intergenic_genes['End'] - intergenic_genes['Start'] + 1
    intergenic_genes = intergenic_genes[intergenic_genes['Length'] >= min_length]
    
    # Filter 2: Remove intergenic genes that overlap with original CDS genes ONLY
    genes_to_keep = []
    
    for _, intergenic_gene in intergenic_genes.iterrows():
        int_start = intergenic_gene['Start']
        int_end = intergenic_gene['End']
        int_name = intergenic_gene['Gene Name']
        
        overlap_found = False
        
        # Check overlap with all original genes
        for _, original_gene in original_genes.iterrows():
            orig_start = original_gene['Start']
            orig_end = original_gene['End']
            orig_name = original_gene['Gene Name']
            orig_category = original_gene.get('Category', '')
            orig_feature_type = original_gene.get('feature type', '')
            
            # Only check overlaps with CDS features
            if orig_feature_type != 'CDS':
                continue  # Skip non-CDS features (origins, ncRNAs, etc.)
            
            # Calculate overlap
            overlap_start = max(int_start, orig_start)
            overlap_end = min(int_end, orig_end)
            overlap_length = max(0, overlap_end - overlap_start + 1)
            
            # Check if significant overlap exists
            if overlap_length >= overlap_threshold:
                # Exception for mobile elements
                if orig_category == 'Mobile Element':
                    continue  # Skip this overlap, don't mark as overlap_found
                else:
                    overlap_found = True
                    break
            
            # Also check if same gene name in close proximity (but skip mobile elements)
            if int_name == orig_name and orig_category != 'Mobile Element':
                distance = min(abs(int_start - orig_end), abs(orig_start - int_end))
                if distance < 500:  # Within 500bp
                    overlap_found = True
                    break
        
        if not overlap_found:
            genes_to_keep.append(intergenic_gene)
    
    # Create filtered intergenic dataframe
    if genes_to_keep:
        filtered_intergenic = pd.DataFrame(genes_to_keep)
        # Remove the temporary Length column
        if 'Length' in filtered_intergenic.columns:
            filtered_intergenic = filtered_intergenic.drop('Length', axis=1)
    else:
        filtered_intergenic = pd.DataFrame()
    
    # Combine filtered results
    if not filtered_intergenic.empty:
        final_df = pd.concat([original_genes, filtered_intergenic], ignore_index=True)
        final_df = final_df.sort_values('Start').reset_index(drop=True)
    else:
        final_df = original_genes.sort_values('Start').reset_index(drop=True)
    
    return final_df


def detect_intergenic_genes_simple(final_dataframe_fixed, sequence, blast_temp_dir):
    """Main function - CLEAN SINGLE BLAST VERSION - SIMPLIFIED OUTPUT"""
    print("🔍 Analyzing intergenic regions to find out additional genes")
    
    try:
        blast_db_prefix = "database_blast/translations_db"
        database_csv_path = "Database/Database.csv"
        
        # Count original genes (excluding non-CDS features)
        original_cds_count = len(final_dataframe_fixed[final_dataframe_fixed.get('feature type', '') == 'CDS'])
        
        # Step 1: Find intergenic regions
        intergenic_regions = find_intergenic_regions_simple(
            final_dataframe_fixed, sequence, min_gap_size=100
        )
        
        if not intergenic_regions:
            print(f"   Original genes: {original_cds_count}")
            print(f"   Additional genes from intergenic regions: 0")
            print(f"   Final gene count: {original_cds_count}")
            return final_dataframe_fixed
        
        # Step 2: BLAST intergenic regions (single BLAST per region) - SILENT
        blast_intergenic_regions_single(intergenic_regions, blast_db_prefix, blast_temp_dir)
        
        # Step 3: Parse results - SILENT
        new_annotations = parse_intergenic_blast_single(
            blast_temp_dir, database_csv_path, intergenic_regions, sequence
        )
        
        if new_annotations.empty:
            print(f"   Original genes: {original_cds_count}")
            print(f"   Additional genes from intergenic regions: 0")
            print(f"   Final gene count: {original_cds_count}")
            return final_dataframe_fixed
        
        # Step 4: Add flag column to existing annotations
        if 'Intergenic_Flag' not in final_dataframe_fixed.columns:
            final_dataframe_fixed['Intergenic_Flag'] = 0
        
        # Step 5: Combine dataframes
        for col in final_dataframe_fixed.columns:
            if col not in new_annotations.columns:
                new_annotations[col] = ''
        
        for col in new_annotations.columns:
            if col not in final_dataframe_fixed.columns:
                final_dataframe_fixed[col] = ''
        
        combined_df = pd.concat([final_dataframe_fixed, new_annotations], ignore_index=True)
        combined_df = combined_df.sort_values('Start').reset_index(drop=True)
        
        # Step 6: Filter overlapping and short genes - SILENT
        filtered_df = filter_overlapping_intergenic_genes(
            combined_df, 
            min_length=90,        # Remove genes < 90bp
            overlap_threshold=50  # Remove genes with >50bp overlap
        )
        
        final_intergenic_count = len(filtered_df[filtered_df['Intergenic_Flag'] == 1])
        final_total_cds = len(filtered_df[filtered_df.get('feature type', '') == 'CDS'])
        
        print(f"   Original genes: {original_cds_count}")
        print(f"   Additional genes from intergenic regions: {final_intergenic_count}")
        print(f"   Final gene count: {final_total_cds}")
        
        return filtered_df
        
    except Exception as e:
        print(f"❌ Error in intergenic detection: {e}")
        return final_dataframe_fixed