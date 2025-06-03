import subprocess
import pandas as pd
from Bio import SeqIO
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from pathlib import Path
from collections import Counter
from Bio import SeqIO, Entrez

import pandas as pd
import subprocess
import re
import os
from Bio import SeqIO

def CDS_fasta_using_prodigal(fasta_path, output_name):
    """
    Runs Prodigal on a FASTA file, parses CDS entries using regex from GBK,
    and returns a DataFrame with Start, End, Strand, Sequence.
    The temporary GBK file is deleted after processing.
    """
    gbk_path = fasta_path + ".gbk"

    # Step 1: Read original sequence from FASTA
    try:
        record = SeqIO.read(fasta_path, "fasta")
        full_seq = record.seq
    except Exception as e:
        print(f"❌ Error reading FASTA file: {e}")
        return

    # Step 2: Run Prodigal
    try:
        subprocess.run([
            "prodigal",
            "-i", fasta_path,
            "-o", gbk_path,
            "-f", "gbk",
            "-q",
            "-p", "meta",
            "-c"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print("❌ Prodigal failed:", e)
        return

    # Step 3: Parse CDS entries from GBK file using regex
    with open(gbk_path, "r") as f:
        gbk_content = f.read()

    cds_pattern = re.compile(r'^\s+CDS\s+(complement\()?(\d+)\.\.(\d+)', re.MULTILINE)
    cds_entries = []

    for match in cds_pattern.finditer(gbk_content):
        is_complement = match.group(1) is not None
        start = int(match.group(2)) - 1  # 0-based index
        end = int(match.group(3))       # end is exclusive for slicing
        strand = '-' if is_complement else '+'

        subseq = full_seq[start:end]
        if strand == '-':
            subseq = subseq.reverse_complement()

        cds_entries.append({
            "Start": start + 1,  # 1-based
            "End": end,
            "Strand": strand,
            "Sequence": str(subseq)
        })

    # Step 4: Delete GBK file
    try:
        os.remove(gbk_path)
    except Exception as e:
        print(f"⚠️ Warning: Failed to delete temporary GBK file: {e}")

    # Step 5: Output
    if cds_entries:
        df = pd.DataFrame(cds_entries)
        print(f"✅ Extracted {len(df)} CDS entries from {output_name} plasmid")
    else:
        print("⚠️ No CDS entries found.")
        df = pd.DataFrame()

    return df





def CDS_genbank_overwrite(
    genbank_file, output_name,
    temp_dir,
    email="hislam2@ur.rochester.edu"
):
    """
    Uses a GenBank file as a FASTA input for Prodigal-based CDS extraction.
    If the GenBank lacks sequence, tries to fetch from NCBI.
    """
    Entrez.email = email
    sequence = None
    accession = None
    temp_fasta = os.path.join(temp_dir, "temp_for_prodigal.fasta")

    # Step 1: Try reading the GenBank file
    try:
        record = SeqIO.read(genbank_file, "genbank")
        accession = record.id
        if record.seq and "N" not in str(record.seq):
            sequence = record.seq
            print(f"📦 Sequence found in GenBank file for accession {accession}")
    except Exception as e:
        print(f"⚠️ Could not read GenBank file: {e}")

    # Step 2: Fallback to NCBI
    if not sequence:
        print("🔍 Sequence not found in GenBank. Attempting to fetch from NCBI...")
        if not accession:
            with open(genbank_file) as f:
                content = f.read()
            match = re.search(r"ACCESSION\s+([A-Z_]+\d+)", content)
            if match:
                accession = match.group(1)
            else:
                print("❌ Could not identify accession.")
                return None, None, None

        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            fasta_record = SeqIO.read(handle, "fasta")
            handle.close()
            sequence = fasta_record.seq
            print("✅ Sequence successfully fetched from NCBI.")
        except Exception as e:
            print("❌ Could not fetch FASTA sequence from NCBI.")
            return None, None, None

    # Step 3: Write sequence to temporary FASTA
    try:
        with open(temp_fasta, "w") as f:
            f.write(f">{accession}\n{sequence}\n")
        print(f"📝 Temporary FASTA written to {temp_fasta}")
    except Exception as e:
        print(f"❌ Failed to write temporary FASTA: {e}")
        return None, None, None

    # Step 4: Run Prodigal
    try:
        df = CDS_fasta_using_prodigal(temp_fasta, output_name)
    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

    return df, str(sequence), len(sequence)


def generate_orf_annotation(cds_df, blast_folder):
    fields = ["Gene Name", "Product", "Category"]

    codon_table_11 = {
        'ttt':'F','ttc':'F','tta':'L','ttg':'L','ctt':'L','ctc':'L','cta':'L','ctg':'L',
        'att':'I','atc':'I','ata':'I','atg':'M',
        'gtt':'V','gtc':'V','gta':'V','gtg':'V',
        'tct':'S','tcc':'S','tca':'S','tcg':'S','agt':'S','agc':'S',
        'cct':'P','ccc':'P','cca':'P','ccg':'P',
        'act':'T','acc':'T','aca':'T','acg':'T',
        'gct':'A','gcc':'A','gca':'A','gcg':'A',
        'tat':'Y','tac':'Y',
        'cat':'H','cac':'H',
        'caa':'Q','cag':'Q',
        'aat':'N','aac':'N',
        'aaa':'K','aag':'K',
        'gat':'D','gac':'D',
        'gaa':'E','gag':'E',
        'tgt':'C','tgc':'C',
        'tgg':'W',
        'cgt':'R','cgc':'R','cga':'R','cgg':'R','aga':'R','agg':'R',
        'ggt':'G','ggc':'G','gga':'G','ggg':'G',
        'taa':'*','tag':'*','tga':'*'
    }

    def translate_dna(seq):
        seq = seq.lower().replace("\n", "").replace(" ", "")
        protein = ''
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa = codon_table_11.get(codon, 'X')
            protein += aa
        return protein

    def get_majority_values(df):
        result = {
            "Gene Name": "ORF",
            "Product": "Open reading frame",
            "Category": "Uncharacterized"
        }
        for field in fields:
            if field in df.columns:
                counts = Counter(df[field].dropna())
                if counts:
                    result[field] = counts.most_common(1)[0][0]
        return result

    blast_folder = Path(blast_folder)
    final_records = []

    for idx in range(len(cds_df)):
        blast_file = blast_folder / f"blast_result_{idx}.csv"

        if blast_file.exists() and blast_file.stat().st_size > 0:
            df = pd.read_csv(blast_file)

            for threshold in [80, 60, 40]:
                filtered = df[df["identity"] > threshold]
                if not filtered.empty:
                    values = get_majority_values(filtered)
                    break
            else:
                values = get_majority_values(pd.DataFrame())
        else:
            values = get_majority_values(pd.DataFrame())

        cds_row = cds_df.iloc[idx]
        protein_seq = translate_dna(cds_row["Sequence"])

        final_records.append({
            "Gene Name": values["Gene Name"],
            "Product": values["Product"],
            "Start": cds_row["Start"],
            "End": cds_row["End"],
            "Strand": cds_row["Strand"],
            "Category": values["Category"],
            "Translation": protein_seq
        })

    final_df = pd.DataFrame(final_records)
    #final_df.to_csv(output_file, index=False)
    #print(f"✅ Final annotated ORFs written to {output_file}")
    return final_df


def combine_features(cds_df, ncrna_df, full_fasta_str):
    # Parse full genome sequence from FASTA
    def parse_full_fasta(fasta_str):
        if fasta_str is None:  # ✨ Handle None sequence
            return ""
        seq_lines = []
        for line in fasta_str.strip().splitlines():
            if not line.startswith(">"):
                seq_lines.append(line.strip())
        return "".join(seq_lines).upper()

    def get_subsequence(seq, start, end, strand):
        try:
            if not seq:  # ✨ Handle empty sequence
                return ""
            subseq = seq[int(start)-1:int(end)]  # 1-based inclusive
            if strand == "-":
                complement = str.maketrans("ACGT", "TGCA")
                subseq = subseq.translate(complement)[::-1]
            return subseq
        except Exception as e:
            print(f"⚠️ Error extracting subsequence: {e}")
            return ""

    genome_seq = parse_full_fasta(full_fasta_str)

    # Process CDS
    cds_df = cds_df.copy()
    cds_df["feature type"] = "CDS"

    # Process ncRNA with extracted sequences
    # Check if ncRNA dataframe is empty
    if ncrna_df.empty:
        print("📝 No ncRNA features found, skipping ncRNA processing")
        # Create empty ncRNA dataframe with correct structure
        ncrna_df = pd.DataFrame(columns=["Gene Name", "Product", "Start", "End", "Strand", "Category", "Translation", "feature type"])
    else:
        ncrna_df = ncrna_df.copy()
        
        # Rename columns safely
        column_mapping = {
            "ncRNA Family": "Gene Name",
            "Description": "Product"
        }
        
        # Only rename columns that exist
        for old_col, new_col in column_mapping.items():
            if old_col in ncrna_df.columns:
                ncrna_df = ncrna_df.rename(columns={old_col: new_col})
        
        ncrna_df["feature type"] = "NC_RNA"
        ncrna_df["Category"] = "Non coding RNA/Regulatory elements"
        
        # Add Translation column with proper error handling
        try:
            # Use a safer approach for adding the Translation column
            translation_sequences = []
            for idx, row in ncrna_df.iterrows():
                seq = get_subsequence(genome_seq, row["Start"], row["End"], row["Strand"])
                translation_sequences.append(seq)
            
            ncrna_df["Translation"] = translation_sequences
            
        except Exception as e:
            print(f"❌ Error adding Translation to ncRNA: {e}")
            print(f"🔍 ncRNA row causing issue: {row}")
            # Add empty Translation column as fallback
            ncrna_df["Translation"] = ""

    # Harmonize columns
    columns_order = ["Gene Name", "Product", "Start", "End", "Strand", "Category", "Translation", "feature type"]
    
    # Ensure all required columns exist in both dataframes
    for col in columns_order:
        if col not in cds_df.columns:
            cds_df[col] = ""
        if col not in ncrna_df.columns:
            ncrna_df[col] = ""
    
    # Select columns in correct order
    cds_df = cds_df[columns_order]
    ncrna_df = ncrna_df[columns_order]

    # Combine all
    combined_df = pd.concat([cds_df, ncrna_df], ignore_index=True)
    
    return combined_df


import pandas as pd
from Bio.Seq import Seq

def Fixing_dataframe(df, fasta_sequence, output_path=None):
    # Load CSV and sanitize column names
    df.columns = df.columns.str.strip()

    # ✨ Handle None sequence
    if fasta_sequence is None:
        print("⚠️ No genomic sequence available for coordinate validation")
        genome_seq = None
    else:
        genome_seq = Seq(fasta_sequence)

    # Validate required columns
    required_cols = ['Category', 'Start', 'End', 'Strand', 'Gene Name']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: '{col}'")

    # Step 1: Merge fragmented mobile elements using 'Gene Name' as identifier
    # ✨ Only if we have sequence (mobile elements need coordinates validation)
    if genome_seq is not None:
        merged_rows = []
        skip_indices = set()

        for i in range(len(df) - 1):
            if i in skip_indices:
                continue

            row = df.iloc[i]
            next_row = df.iloc[i + 1]

            if (row['Category'] == 'Mobile Element' and
                next_row['Category'] == 'Mobile Element' and
                row['Gene Name'] == next_row['Gene Name'] and
                next_row['Start'] - row['End'] <= 1000 and
                row['Strand'] == next_row['Strand']):

                new_start = row['Start']
                new_end = next_row['End']
                strand = row['Strand']
                merged_seq = genome_seq[new_start - 1:new_end]
                if strand == "-":
                    merged_seq = merged_seq.reverse_complement()

                new_row = row.copy()
                new_row['Start'] = new_start
                new_row['End'] = new_end
                new_row['Translation'] = str(merged_seq)
                merged_rows.append(new_row)
                skip_indices.update([i, i + 1])
            else:
                if i not in skip_indices:
                    merged_rows.append(row)

        if len(df) - 1 not in skip_indices:
            merged_rows.append(df.iloc[-1])

        df = pd.DataFrame(merged_rows).reset_index(drop=True)

    # Step 2: Apply corrections
    for idx, row in df.iterrows():
        gene = str(row.get("Gene Name", "")).strip()
        prod = str(row.get("Product", "")).strip()
        ftype = str(row.get("feature type", "")).strip()
        cat = str(row.get("Category", "")).strip()

        # oriT
        if gene.lower() == "orit":
            df.at[idx, "Product"] = "Origin of Transfer"
            df.at[idx, "Category"] = "Origin of Transfer"
            df.at[idx, "feature type"] = "ORIT"
            
        # oriV
        elif gene.lower() == "oriv":
            df.at[idx, "Product"] = "Origin of Replication"
            df.at[idx, "Category"] = "Origin of Replication"
            df.at[idx, "feature type"] = "ORIV"

        # Mobile Element → MGE
        if cat == "Mobile Element":
            df.at[idx, "feature type"] = "MGE"

        if cat =="Replication":
            df.at[idx, "feature type"] = "Replicon"
            df.at[idx, "Category"] = "Replicon"
            df.at[idx, "Product"] = "Predicted replicon"

        if cat == "Uncharacterized":
            df.at[idx, "Category"] = "Open reading frame"

        if ftype=="NC_RNA":
            df.at[idx,"Category"]="Non coding RNA/Regulatory elements"

        # Fix NC_RNA Translation only if we have genomic sequence
        if ftype == "NC_RNA" and genome_seq is not None:
            try:
                s = int(row["Start"])
                e = int(row["End"])
                start, end = min(s, e) - 1, max(s, e)
                strand = str(row["Strand"]).strip()
                seq = genome_seq[start:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                df.at[idx, "Translation"] = str(seq)
            except Exception as e:
                df.at[idx, "Translation"] = f"[Error: {e}]"

    for idx, row in df.iterrows():
        if row.get("Category") == "Replicon":
            gene_name = str(row.get("Gene Name", ""))
            if ">" in gene_name:
                gene_core = gene_name.split(">")[1]
                main_part = gene_core.split("_")[0]
                df.at[idx, "Gene Name"] = main_part

    if output_path:
        df.to_csv(output_path, index=False)

    return df


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import pandas as pd
from datetime import datetime

def write_genbank_from_annotation(df, sequence, output_path, record_id, record_name, description):
    """
    Converts an annotation DataFrame and nucleotide sequence string into a GenBank file.
    ✨ Handles None sequence by creating dummy sequence based on max coordinate
    """
    
    # ✨ Handle None sequence
    if sequence is None:
        # Calculate required sequence length from max end coordinate
        max_position = int(df['End'].max()) if not df.empty else 1000
        sequence = 'N' * max_position  # Create dummy sequence of N's
        print(f"⚠️ Using placeholder sequence ({max_position} bp) for GenBank generation")
    
    seq = Seq(sequence)
    
    record = SeqRecord(
        seq,
        id=record_id,
        name=record_name,
        description=description,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",  # or "linear"
            "date": datetime.today().strftime("%d-%b-%Y").upper(),
            "data_file_division": "BCT"  # For bacterial plasmids
        }
    )

    # Add source feature
    source_feature = SeqFeature(
        FeatureLocation(0, len(seq)),
        type="source",
        qualifiers={
            "organism": "synthetic construct",
            "mol_type": "plasmid DNA"
        }
    )
    record.features.append(source_feature)

    for _, row in df.iterrows():
        try:
            start = int(row["Start"])
            end = int(row["End"])
            strand = 1 if row["Strand"].strip() == "+" else -1
            feature_type = str(row["feature type"]).strip()

            # Ensure start <= end
            if start > end:
                start, end = end, start

            # Convert to 0-based indexing for Biopython
            start_0b = start - 1
            end_0b = end
            location = FeatureLocation(start_0b, end_0b, strand=strand)

            qualifiers = {
                "gene": str(row["Gene Name"]),
                "product": str(row["Product"]),
                "category": str(row["Category"]),
            }

            if feature_type == "CDS" and pd.notnull(row.get("Translation", None)):
                qualifiers["translation"] = str(row["Translation"])
            else:
                # ✨ Handle sequence extraction safely
                if start_0b < len(sequence) and end_0b <= len(sequence):
                    dna_seq = sequence[start_0b:end_0b]
                    if strand == -1:
                        dna_seq = str(Seq(dna_seq).reverse_complement())
                    qualifiers["sequence"] = dna_seq
                else:
                    # Coordinates exceed sequence length (shouldn't happen with dummy sequence)
                    qualifiers["sequence"] = "N" * (end_0b - start_0b)

            feature = SeqFeature(
                location=location,
                type=feature_type,
                qualifiers=qualifiers
            )
            record.features.append(feature)

        except Exception as e:
            print(f"⚠️ Skipping row due to error: {e}\nRow:\n{row}\n")

    # Write GenBank file
    SeqIO.write(record, output_path, "genbank")

'''def generate_plasmid_summary_stats(final_dataframe, plasmid_sequence, plasmid_name, plasmid_length):
    """Generate comprehensive summary statistics for a plasmid"""
    
    print(f"\n📊 {plasmid_name} - Summary Statistics")
    print("=" * 50)
    
    # Basic sequence statistics
    if plasmid_sequence:
        gc_count = plasmid_sequence.upper().count('G') + plasmid_sequence.upper().count('C')
        gc_content = (gc_count / len(plasmid_sequence)) * 100
        at_content = 100 - gc_content
    else:
        gc_content = "N/A"
        at_content = "N/A"
    
    print(f"🧬 SEQUENCE STATISTICS:")
    print(f"   Length: {plasmid_length:,} bp")
    print(f"   GC content: {gc_content:.1f}%" if isinstance(gc_content, float) else f"   GC content: {gc_content}")
    print(f"   AT content: {at_content:.1f}%" if isinstance(at_content, float) else f"   AT content: {at_content}")
    
    # Gene content statistics
    total_genes = len(final_dataframe)
    cds_genes = len(final_dataframe[final_dataframe['feature type'] == 'CDS'])
    
    # Count annotated vs hypothetical
    annotations = final_dataframe.get('Product', final_dataframe.get('Annotation', []))
    hypothetical_count = sum(1 for ann in annotations if 'hypothetical' in str(ann).lower())
    annotated_count = cds_genes - hypothetical_count
    
    # Calculate coding density
    if plasmid_sequence and cds_genes > 0:
        total_coding_length = 0
        for _, row in final_dataframe.iterrows():
            if row['feature type'] == 'CDS':
                gene_length = abs(row['End'] - row['Start']) + 1
                total_coding_length += gene_length
        coding_density = (total_coding_length / plasmid_length) * 100
        avg_gene_length = total_coding_length / cds_genes
        gene_density = (cds_genes / plasmid_length) * 1000  # genes per kb
    else:
        coding_density = "N/A"
        avg_gene_length = "N/A"
        gene_density = "N/A"
    
    print(f"\n🧬 GENE CONTENT:")
    print(f"   Total features: {total_genes}")
    print(f"   CDS genes: {cds_genes}")
    print(f"   Annotated genes: {annotated_count} ({annotated_count/cds_genes*100:.1f}%)" if cds_genes > 0 else "   Annotated genes: 0")
    print(f"   Hypothetical proteins: {hypothetical_count} ({hypothetical_count/cds_genes*100:.1f}%)" if cds_genes > 0 else "   Hypothetical proteins: 0")
    print(f"   Coding density: {coding_density:.1f}%" if isinstance(coding_density, float) else f"   Coding density: {coding_density}")
    print(f"   Average gene length: {avg_gene_length:.0f} bp" if isinstance(avg_gene_length, float) else f"   Average gene length: {avg_gene_length}")
    print(f"   Gene density: {gene_density:.1f} genes/kb" if isinstance(gene_density, float) else f"   Gene density: {gene_density}")
    
    # Functional category counts
    functional_categories = {
        'AMR': 0,
        'Virulence': 0, 
        'Mobilization': 0,
        'Replication': 0,
        'Metabolism': 0,
        'Regulatory': 0,
        'Transposon': 0
    }
    
    # Count genes by function based on annotations
    for _, row in final_dataframe.iterrows():
        if row['feature type'] != 'CDS':
            continue
            
        annotation = str(row.get('Product', row.get('Annotation', ''))).lower()
        category = str(row.get('Category', '')).lower()
        
        # AMR keywords
        if any(keyword in annotation for keyword in ['resistance', 'beta-lactam', 'antibiotic', 'efflux', 'bla', 'tet', 'sul', 'qnr']):
            functional_categories['AMR'] += 1
        # Virulence keywords
        elif any(keyword in annotation for keyword in ['virulence', 'toxin', 'adhesin', 'hemolysin']):
            functional_categories['Virulence'] += 1
        # Mobilization keywords
        elif any(keyword in annotation for keyword in ['conjugation', 'transfer', 'tra', 'trb', 'mobilization']):
            functional_categories['Mobilization'] += 1
        # Replication keywords
        elif any(keyword in annotation for keyword in ['replication', 'rep', 'origin']) or 'replication' in category:
            functional_categories['Replication'] += 1
        # Metabolism keywords
        elif any(keyword in annotation for keyword in ['metabolism', 'transport', 'kinase', 'synthase']):
            functional_categories['Metabolism'] += 1
        # Regulatory keywords
        elif any(keyword in annotation for keyword in ['regulation', 'transcription', 'regulator', 'repressor']):
            functional_categories['Regulatory'] += 1
        # Transposon keywords
        elif any(keyword in annotation for keyword in ['transposon', 'insertion', 'transposase']):
            functional_categories['Transposon'] += 1
    
    print(f"\n🎯 FUNCTIONAL CATEGORIES:")
    for category, count in functional_categories.items():
        percentage = (count / cds_genes * 100) if cds_genes > 0 else 0
        print(f"   {category}: {count} ({percentage:.1f}%)")
    
    # Mobile genetic elements (count from feature type)
    mobile_elements = {
        'oriT': len(final_dataframe[final_dataframe['feature type'] == 'oriT']),
        'oriC': len(final_dataframe[final_dataframe['feature type'] == 'oriC']),
        'Transposon': len(final_dataframe[final_dataframe['feature type'] == 'transposon']),
        'Replicon': len(final_dataframe[final_dataframe['feature type'] == 'replicon']),
        'ncRNA': len(final_dataframe[final_dataframe['feature type'] == 'ncRNA'])
    }
    
    print(f"\n🚀 MOBILE GENETIC ELEMENTS:")
    for element, count in mobile_elements.items():
        print(f"   {element}: {count}")
    
    # Calculate some ratios
    print(f"\n📈 KEY RATIOS:")
    if cds_genes > 0:
        amr_ratio = functional_categories['AMR'] / cds_genes * 100
        mobile_ratio = functional_categories['Mobilization'] / cds_genes * 100
        known_ratio = annotated_count / cds_genes * 100
        print(f"   AMR gene ratio: {amr_ratio:.1f}%")
        print(f"   Mobilization gene ratio: {mobile_ratio:.1f}%") 
        print(f"   Known function ratio: {known_ratio:.1f}%")
    
    print("=" * 50)
    
    # Return summary dict for potential CSV export
    summary_dict = {
        'Plasmid_Name': plasmid_name,
        'Length_bp': plasmid_length,
        'GC_Content': gc_content if isinstance(gc_content, float) else None,
        'Total_Genes': total_genes,
        'CDS_Genes': cds_genes,
        'Annotated_Genes': annotated_count,
        'Hypothetical_Genes': hypothetical_count,
        'Coding_Density': coding_density if isinstance(coding_density, float) else None,
        'Average_Gene_Length': avg_gene_length if isinstance(avg_gene_length, float) else None,
        'Gene_Density_per_kb': gene_density if isinstance(gene_density, float) else None,
        'AMR_Genes': functional_categories['AMR'],
        'Virulence_Genes': functional_categories['Virulence'],
        'Mobilization_Genes': functional_categories['Mobilization'],
        'Replication_Genes': functional_categories['Replication'],
        'oriT_Count': mobile_elements['oriT'],
        'Transposon_Count': mobile_elements['Transposon'],
        'ncRNA_Count': mobile_elements['ncRNA']
    }
    
    return summary_dict


def save_summary_stats_csv(summary_dict, output_folder, plasmid_name):
    """Save summary statistics to CSV file"""
    import pandas as pd
    
    summary_df = pd.DataFrame([summary_dict])
    csv_file = os.path.join(output_folder, f"{plasmid_name}_summary_stats.csv")
    summary_df.to_csv(csv_file, index=False)
    print(f"📊 Summary statistics saved: {csv_file}")
    return csv_file'''

def print_sequence_statistics(plasmid_sequence, plasmid_name, plasmid_length):
    """Print basic sequence statistics for a plasmid"""
    
    # Only print if sequence is available
    if plasmid_sequence:
        print(f"\n🧬 SEQUENCE STATISTICS:")
        print(f"   Length: {plasmid_length:,} bp")
        
        # Calculate GC content
        sequence_upper = plasmid_sequence.upper()
        gc_count = sequence_upper.count('G') + sequence_upper.count('C')
        gc_content = (gc_count / len(plasmid_sequence)) * 100
        at_content = 100 - gc_content
        
        print(f"   GC content: {gc_content:.1f}%")
        print(f"   AT content: {at_content:.1f}%")