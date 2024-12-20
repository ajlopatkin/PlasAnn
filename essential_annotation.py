import os
import subprocess as sp
from Bio.Seq import Seq
import pandas as pd
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from io import StringIO
import re


def getsequencelength(plasmid,pathofdir):
    length_of_fasta= sp.getoutput(f"awk '{{if(NR>1)print}}' {pathofdir}/{plasmid}.fasta | awk '{{print length}}' | awk '{{n += $1}}; END{{print n}}'")
    return length_of_fasta

def getthesequences(plasmid,pathofdir):
    only_sequence = sp.getoutput(f"awk '{{if(NR>1)print}}' {pathofdir}/{plasmid}.fasta | tr -d '\r\n' ")
    return only_sequence

def getpositionsofCDS(plasmid):
    os.system(f"grep CDS tmp_files/{plasmid}prodigal.gbk | grep -o '[[:digit:]]*' > tmp_files/positions.txt")  #this one gets all the positions of CDS
    with open("tmp_files/positions.txt", "r") as file_obj:
        file_data = file_obj.read() 
        lines = file_data.splitlines() 
        list_of_positions = [(lines[i], lines[i+1]) for i in range(0, len(lines)-1,2)]
    #os.remove("../tmp_files/positions.txt")
    return list_of_positions


# This one gets the complement strat and end positions 
def complementpositions(plasmid):
    os.system(f"grep complement  tmp_files/{plasmid}prodigal.gbk |awk '{{print $2}}'|grep -o '([^)]*)'| awk -F\\. '{{print $1}}'| grep -o '[[:digit:]]*' > tmp_files/complementstarting.txt")  #This one gets all the complement positions [starting]
    with open("tmp_files/complementstarting.txt","r") as file:
        list_of_complementstart = []
        for line in file:
            line = line.strip()
            list_of_complementstart.append(line)
    #os.remove("../tmp_files/complementstarting.txt")

    os.system(f"grep complement  tmp_files/{plasmid}prodigal.gbk |awk '{{print $2}}'|grep -o '([^)]*)'| awk -F\\. '{{print $3}}'| grep -o '[[:digit:]]*' > tmp_files/complementending.txt")   #This one gets all the complement positions [ending]
    with open("tmp_files/complementending.txt","r") as file:
        list_of_complementend = []
        for line in file:
            line = line.strip()
            list_of_complementend.append(line)
    #os.remove("../tmp_files/complementending.txt")
    return list_of_complementstart, list_of_complementend


def getlistofDNACDS(list_of_positions,full_sequence):
    list_of_cds = []
    for i in range(len(list_of_positions)):
        x = int(list_of_positions[i][0])-1
        y = int(list_of_positions[i][1])-1
        list_of_cds.append(full_sequence[x:y])
    return list_of_cds


#The one that makes a query file from the list of coding sequences
def makequeryfastafordbsearch(list_of_cds):
    query_file = open(r'tmp_files/Query_Fasta.fsa', 'w+')
    out = '\n'.join(['>Coding Sequence' + str(i+1) + "\n" + j for i,j in enumerate(list_of_cds)])
    query_file.write(out)
    query_file.close()

def makedatabasefromcsvfile(database_name):
    database = pd.read_csv(database_name,low_memory=False)
    database['identifier'] = database['Gene Name'] + '#"' + database.index.astype(str) + '"'
    target_file = open(r'tmp_files/target_fasta.fsa','w+')
    out1 = '\n'.join(['>'+str(database['identifier'][i])+"\n"+str(database['Translation'][i]) for i in range(len(database))])
    target_file.write(out1)
    target_file.close()
    if not os.path.exists("makedb_folder/blastdb"):
        os.makedirs("makedb_folder/blastdb")
    os.system('makeblastdb -in tmp_files/target_fasta.fsa -dbtype prot -out makedb_folder/blastdb/custom_database > /dev/null 2>&1')
    return database



def extract_gene_info(gene_string):
    if pd.isna(gene_string):
        return None, None
    try:
        gene_name, index_part = gene_string.split('#"')
        index_part = index_part.rstrip('"')  # Remove the trailing quote
        index = int(index_part)  # Convert index to integer
        return gene_name, index
    except ValueError:
        return None, None


def select_gene_row(df):
    # First try to filter rows with 'Percent Identity' = 100
    filtered_df = df[df['Pident'] >= 85]
    
    # If no rows with 100% identity, or less than 10 rows with 100% identity, consider rows with 'Percent Identity' > 90%
    #if filtered_df.empty or len(filtered_df) < 5:
    if filtered_df.empty:
        filtered_df = df[df['Pident'] > 70]
    
    # Count the occurrences of each gene
    gene_counts = filtered_df['Gene Name'].value_counts()
    max_occurrence = gene_counts.max()
    # Identify genes with the maximum occurrence
    most_frequent_genes = gene_counts[gene_counts == max_occurrence].index.tolist()
    # Filter for the most frequent genes
    most_frequent_df = filtered_df[filtered_df['Gene Name'].isin(most_frequent_genes)]
    # Select the row with the lowest E-value for each gene
    selected_rows = most_frequent_df.loc[most_frequent_df.groupby('Gene Name')['E-value'].idxmin()]
    return selected_rows

def remove_similar_positions(df):
    df['Start Position Range'] = df['Start Position'] // 50
    df['End Position Range'] = df['End Position'] // 50
    
    # Group by Start Position Range and End Position Range and select row with minimum E-value in each group
    df = df.loc[df.groupby(['Start Position Range', 'End Position Range'])['E-value'].idxmin()]

    # Drop the temporary columns used for grouping
    df = df.drop(columns=['Start Position Range', 'End Position Range'])
    
    return df

'''def initial_blast_against_database(list_of_cds,list_of_positions,database):

    final_dataframe = pd.DataFrame()
    name_of_query = []
    titlelist = []
    lengthlist = []
    scorelist = []
    gaplist = []
    evallist = []
    sequencelist = []
    startpos = []
    endpos = []
    pident_value = []  # Ensure this list is initialized
    query_lengths = []  # Ensure this list is initialized
    subject_lengths = []  # Ensure this list is initialized

    # Assuming list_of_cds and list_of_positions are defined somewhere above this snippet
    for i in range(len(list_of_cds)):
        with open(r'../tmp_files/query1.fsa', 'w+') as query1_file:
            out = '\n'.join(['>Coding Sequence' + str(i + 1) + "\n" + list_of_cds[i]])
            query1_file.write(out)

        # Ensure to give the correct path to your blastx executable
        os.system('/Users/habibulislam/ncbi-blast-2.15.0+/bin/blastx -query ../tmp_files/query1.fsa -db ../makedb_folder/blastdb/custom_database -outfmt 5 -out ../tmp_files/result1.xml')
        with open("../tmp_files/result1.xml", "r") as result:
            records = NCBIXML.parse(result)
            item = next(records)
            
            for alignment in item.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.001 and hsp.identities / alignment.length * 100 > 90:
                        name_of_query.append("Coding Sequence" + str(i + 1))
                        startpos.append(list_of_positions[i][0])
                        endpos.append(list_of_positions[i][1])
                        titlelist.append(alignment.title.split(' ')[1])
                        lengthlist.append(alignment.length)
                        scorelist.append(hsp.score)
                        pident_value.append(hsp.identities / alignment.length * 100)
                        gaplist.append(hsp.gaps)
                        evallist.append(hsp.expect)
                        sequencelist.append(hsp.query)
                        query_lengths.append(len(list_of_cds[i]))  # Assuming this is the intended measure for query_length
                        subject_lengths.append(alignment.length)

    # Check if all lists have the same length before creating the DataFrame
    #print({ "List Name": len(lst) for lst in [name_of_query, titlelist, lengthlist, scorelist, gaplist, evallist, sequencelist, startpos, endpos, pident_value, query_lengths, subject_lengths] })

    # Construct the DataFrame from the collected lists
    final_dataframe = pd.DataFrame({
        'Name of Query': name_of_query,
        'Start Position': startpos,
        'End Position': endpos,
        'Title': titlelist,
        'Length': lengthlist,
        'Score': scorelist,
        'Gaps': gaplist,
        'E-value': evallist,
        'Sequence': sequencelist,
        'Pident': pident_value,
        'Query_Length': query_lengths,
        'Subject_Length': subject_lengths
    })
    final_dataframe['Gene Name'], final_dataframe['index'] = zip(*final_dataframe['Title'].apply(extract_gene_info))
    final_dataframe.dropna(subset=['Gene Name', 'index'], inplace=True)
    unique_query_names = final_dataframe['Name of Query'].unique()
    # Apply the selection process to each unique 'QueryName' and combine the results
    selected_rows_per_query = pd.concat([select_gene_row(final_dataframe[final_dataframe['Name of Query'] == query_name]) for query_name in unique_query_names])
    selected_rows_per_query = selected_rows_per_query.reset_index(drop=True)
    Product_list=[]
    for i in range(len(selected_rows_per_query)):
        tmp=int(selected_rows_per_query['index'][i])
        Product_list.append(database['Category'][tmp])
    selected_rows_per_query['Product']=Product_list
    selected_rows_per_query['index']=selected_rows_per_query['index'].astype(int)
    Product_list=[]
    for i in range(len(selected_rows_per_query)):
        tmp_index=int(selected_rows_per_query['index'][i])
        Product_list.append(database['Product'][tmp_index])
        #Ab_flag.append(database['Flag'][tmp_index])
    selected_rows_per_query['Original Product']=Product_list
    selected_rows_per_query= selected_rows_per_query.rename(columns={"Product": "Category"})
    selected_rows_per_query= selected_rows_per_query.rename(columns={"Original Product": "Product"})
    columns = ['index','Title', 'Pident', 'Query_Length','Subject_Length']
    selected_rows_per_query.drop(columns, inplace=True, axis=1)
    return selected_rows_per_query'''

def run_blast(sequence, index, list_of_positions):
    query_path = f'tmp_files/query_{index}.fsa'
    result_path = f'tmp_files/result_{index}.xml'
    
    # Write the sequence to a temporary FASTA file
    with open(query_path, 'w') as query_file:
        query_file.write(f'>Coding Sequence{index}\n{sequence}')
    
    # Run BLAST
    os.system(f'blastx -query {query_path} -db makedb_folder/blastdb/custom_database -outfmt 5 -out {result_path}')
    
    # Parse BLAST result
    with open(result_path, 'r') as result_file:
        records = NCBIXML.parse(result_file)
        item = next(records)
        
        results = []
        for alignment in item.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.001 and hsp.identities / alignment.length * 100 > 60:
                    results.append({
                        "Name of Query": f"Coding Sequence{index}",
                        "Start Position": list_of_positions[index][0],
                        "End Position": list_of_positions[index][1],
                        "Title": alignment.title.split(' ')[1],
                        "Length": alignment.length,
                        "Score": hsp.score,
                        "Pident": hsp.identities / alignment.length * 100,
                        "Gaps": hsp.gaps,
                        "E-value": hsp.expect,
                        "Sequence": hsp.query,
                        "Query Length": len(sequence),
                        "Subject Length": alignment.length
                    })
        return results

def initial_blast_against_database(list_of_cds, list_of_positions, database):
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_blast, seq, idx, list_of_positions) for idx, seq in enumerate(list_of_cds)]
        results = [future.result() for future in futures]

    # Flatten list of lists
    flat_results = [item for sublist in results for item in sublist]
    
    # Create DataFrame from results
    final_dataframe = pd.DataFrame(flat_results)
    if not final_dataframe.empty:
        final_dataframe['Gene Name'], final_dataframe['index'] = zip(*final_dataframe['Title'].apply(extract_gene_info))
        final_dataframe.dropna(subset=['Gene Name', 'index'], inplace=True)
        
        # Further processing
        unique_query_names = final_dataframe['Name of Query'].unique()
        selected_rows_per_query = pd.concat([select_gene_row(final_dataframe[final_dataframe['Name of Query'] == query_name]) for query_name in unique_query_names])
        selected_rows_per_query.reset_index(drop=True, inplace=True)
        
        # Append product and category information
        selected_rows_per_query['Category'] = [database['Category'][int(index)] for index in selected_rows_per_query['index']]
        selected_rows_per_query['Product'] = [database['Product'][int(index)] for index in selected_rows_per_query['index']]
        
        # Clean up DataFrame
        selected_rows_per_query.drop(['index', 'Title', 'Pident', 'Query Length', 'Subject Length'], axis=1, inplace=True)

    return selected_rows_per_query




def blast_against_oric_dataframe(oric_database,plasmid, path):
    if not os.path.exists("makedb_folder/oricdb"):
        os.makedirs("makedb_folder/oricdb")
        # first make oric database from the fasta file we have. 
    os.system(f'makeblastdb -in {oric_database} -dbtype nucl -out makedb_folder/oricdb/oric_database > /dev/null 2>&1')
    # run blast against that oric databse
    os.system(f'blastn -query {path}/{plasmid}.fasta -db makedb_folder/oricdb/oric_database  -outfmt 5 -out tmp_files/resultoric.xml')

    #make a dataframe from the oric balst result
    import pandas as pd
    from Bio.Blast import NCBIXML

    # Initialize the lists to store information
    name_of_query = []
    titlelist = []
    lengthlist = []
    scorelist = []
    gaplist = []
    evallist = []
    sequencelist = []
    startpos = []
    endpos = []
    product_list = []
    index=[]
    ab_flag=[]

    # Open and parse the BLAST XML result file
    with open("tmp_files/resultoric.xml", "r") as result:
        records = NCBIXML.parse(result)
        item = next(records)
        j = 1  # Start the counter from 1

        for alignment in item.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01 and hsp.identities / alignment.length * 100 > 90:
                    name_of_query.append("OriC Sequence" + str(j))
                    startpos.append(hsp.query_start)
                    endpos.append(hsp.query_end)
                    titlelist.append("oriC")
                    lengthlist.append(alignment.length)
                    scorelist.append(hsp.score)
                    gaplist.append(hsp.gaps)
                    evallist.append(hsp.expect)
                    sequencelist.append(hsp.query)
                    product_list.append('origin of replication')
                    index.append(None),
                    


    # Create the DataFrame from the collected data
    oric_dataframe = pd.DataFrame({
        "Query Name": name_of_query,
        "Gene Name": titlelist,
        "Length": lengthlist,
        "Score": scorelist,
        "Gaps": gaplist,
        "E-value": evallist,
        "Sequence": sequencelist,
        "Start Position": startpos,
        "End Position": endpos,
        "Product": product_list,
        "index":index
    })


    oric_dataframe

    tolerance = 50
    # Sort the dataframe by 'Start Position' and 'End Position'
    df_sorted = oric_dataframe.sort_values(by=['Start Position', 'End Position'])
    rows_to_keep = []
    # Iterate through the sorted dataframe
    for i in range(len(df_sorted)):
        if i == 0:
            rows_to_keep.append(i)
        else:
            current_row = df_sorted.iloc[i]
            previous_row = df_sorted.iloc[rows_to_keep[-1]]
            if abs(current_row['Start Position'] - previous_row['Start Position']) > tolerance or \
            abs(current_row['End Position'] - previous_row['End Position']) > tolerance:
                rows_to_keep.append(i)

    # Filter the dataframe to keep only the identified rows
    oric_dataframe = df_sorted.iloc[rows_to_keep]
    oric_dataframe= oric_dataframe.rename(columns={"Query Name": "Name of Query"})
    oric_dataframe['Category']='Origin of Replication'

    return oric_dataframe

def blast_against_orit_dataframe(orit_database, plasmid, path):
    if not os.path.exists("makedb_folder/oritdb"):
        os.makedirs("makedb_folder/oritdb")
        # first make oric database from the fasta file we have. 
    os.system(f'makeblastdb -in {orit_database} -dbtype nucl -out makedb_folder/oricdb/orit_database > /dev/null 2>&1')
    # run blast against that oric databse
    os.system(f'blastn -query {path}/{plasmid}.fasta -db makedb_folder/oricdb/orit_database  -outfmt 5 -out tmp_files/resultorit.xml')

    #make a dataframe from the oric balst result
    import pandas as pd
    from Bio.Blast import NCBIXML

    # Initialize the lists to store information
    name_of_query = []
    titlelist = []
    lengthlist = []
    scorelist = []
    gaplist = []
    evallist = []
    sequencelist = []
    startpos = []
    endpos = []
    product_list = []
    index=[]


    # Open and parse the BLAST XML result file
    with open("tmp_files/resultorit.xml", "r") as result:
        records = NCBIXML.parse(result)
        item = next(records)
        j = 1  # Start the counter from 1

        for alignment in item.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01 and hsp.identities / alignment.length * 100 > 90:
                    name_of_query.append("Orit Sequence" + str(j))
                    startpos.append(hsp.query_start)
                    endpos.append(hsp.query_end)
                    titlelist.append("oriT")
                    lengthlist.append(alignment.length)
                    scorelist.append(hsp.score)
                    gaplist.append(hsp.gaps)
                    evallist.append(hsp.expect)
                    sequencelist.append(hsp.query)
                    product_list.append('origin of transfer')
                    index.append(None)
                    j += 1  # Increment the counter after appending

    # Create the DataFrame from the collected data
    orit_dataframe = pd.DataFrame({
        "Query Name": name_of_query,
        "Gene Name": titlelist,
        "Length": lengthlist,
        "Score": scorelist,
        "Gaps": gaplist,
        "E-value": evallist,
        "Sequence": sequencelist,
        "Start Position": startpos,
        "End Position": endpos,
        "Product": product_list,
        "index":index,
    })


    tolerance = 50
    # Sort the dataframe by 'Start Position' and 'End Position'
    df_sorted = orit_dataframe.sort_values(by=['Start Position', 'End Position'])
    rows_to_keep = []
    # Iterate through the sorted dataframe
    for i in range(len(df_sorted)):
        if i == 0:
            rows_to_keep.append(i)
        else:
            current_row = df_sorted.iloc[i]
            previous_row = df_sorted.iloc[rows_to_keep[-1]]
            if abs(current_row['Start Position'] - previous_row['Start Position']) > tolerance or \
            abs(current_row['End Position'] - previous_row['End Position']) > tolerance:
                rows_to_keep.append(i)

    # Filter the dataframe to keep only the identified rows
    orit_dataframe = df_sorted.iloc[rows_to_keep]
    orit_dataframe= orit_dataframe.rename(columns={"Query Name": "Name of Query"})
    orit_dataframe['Category']='Origin of Transfer'

    return orit_dataframe


def blast_against_transposon_database(transposon_database,plasmid,path,selected_rows_per_query):
    # now do the same for transposon 

    if not os.path.exists("makedb_folder/transposondb"):
        os.makedirs("makedb_folder/transposondb")
        # first make oric database from the fasta file we have. 
    os.system(f'makeblastdb -in {transposon_database} -dbtype nucl -out makedb_folder/transposondb/transposon_database > /dev/null 2>&1')
    # run blast against that oric databse
    os.system(f'blastn -query {path}/{plasmid}.fasta -db makedb_folder/transposondb/transposon_database  -outfmt 5 -out tmp_files/resulttransposon.xml')

    #make a dataframe from the oric balst result
    import pandas as pd
    from Bio.Blast import NCBIXML

    # Initialize the lists to store information
    name_of_query = []
    titlelist = []
    lengthlist = []
    scorelist = []
    gaplist = []
    evallist = []
    sequencelist = []
    startpos = []
    endpos = []
    product_list = []


    # Open and parse the BLAST XML result file
    with open("tmp_files/resulttransposon.xml", "r") as result:
        records = NCBIXML.parse(result)
        item = next(records)
        j = 1  # Start the counter from 1

        for alignment in item.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01 and hsp.identities / alignment.length * 100 > 90:
                    name_of_query.append("transposon Sequence"+str(j))
                    startpos.append(hsp.query_start)
                    endpos.append(hsp.query_end)
                    titlelist.append(alignment.title.split(' ')[1]) 
                    lengthlist.append(alignment.length)
                    scorelist.append(hsp.score)
                    gaplist.append(hsp.gaps)
                    evallist.append(hsp.expect)
                    sequencelist.append(hsp.query)
                    product_list.append('Mobile Genetic Element')
                    j += 1  # Increment the counter after appending

    # Create the DataFrame from the collected data
    transposon_dataframe = pd.DataFrame({
        "Query Name": name_of_query,
        "Gene Name": titlelist,
        "Length": lengthlist,
        "Score": scorelist,
        "Gaps": gaplist,
        "E-value": evallist,
        "Sequence": sequencelist,
        "Start Position": startpos,
        "End Position": endpos,
        "Product": product_list,
    })

    # Ensure columns exist in the dataframes
    required_columns = ['Start Position', 'End Position', 'Product']
    if not all(col in selected_rows_per_query.columns for col in required_columns):
        raise ValueError(f"One or more required columns are missing in final_dataframe: {required_columns}")

    transposon_required_columns = ['Start Position', 'End Position', 'Length']
    if not all(col in transposon_dataframe.columns for col in transposon_required_columns):
        raise ValueError(f"One or more required columns are missing in transposon_dataframe: {transposon_required_columns}")

    # Convert 'Start Position' and 'End Position' columns to integers
    annotation_df = selected_rows_per_query.copy()
    transposon_df = transposon_dataframe.copy()

    annotation_df['Start Position'] = annotation_df['Start Position'].astype(int)
    annotation_df['End Position'] = annotation_df['End Position'].astype(int)
    transposon_df['Start Position'] = transposon_df['Start Position'].astype(int)
    transposon_df['End Position'] = transposon_df['End Position'].astype(int)

    # Filter the annotation_df for rows where Abricate Flag is 1
    annotation_df_filtered = annotation_df[annotation_df['Product'] == 'Antibiotic resistance']

    # Initialize an empty dataframe to store the matching rows from transposon_df
    matching_transposon_df = pd.DataFrame()

    if not annotation_df_filtered.empty:
        # Iterate over the filtered annotation_df to find matching rows in transposon_df
        for _, row in annotation_df_filtered.iterrows():
            start_pos = row['Start Position']
            end_pos = row['End Position']
            matching_rows = transposon_df[(transposon_df['Start Position'] <= start_pos) & (transposon_df['End Position'] >= end_pos)]
            matching_transposon_df = pd.concat([matching_transposon_df, matching_rows], ignore_index=True)

        # Display the matching rows
        if not matching_transposon_df.empty:
            unique_transposon_df = matching_transposon_df.loc[matching_transposon_df.groupby('Start Position')['Length'].idxmax()]

            # Function to check for overlaps and remove rows with lower length in case of overlap
            def remove_overlapping_rows(df):
                df_sorted = df.sort_values(by=['Start Position']).reset_index(drop=True)
                non_overlapping_rows = []

                current_row = df_sorted.iloc[0]
                non_overlapping_rows.append(current_row)

                for index in range(1, len(df_sorted)):
                    next_row = df_sorted.iloc[index]
                    if next_row['Start Position'] > current_row['End Position']:
                        non_overlapping_rows.append(next_row)
                        current_row = next_row

                return pd.DataFrame(non_overlapping_rows)

            # Apply the function to the dataframe
            non_overlapping_transposon_df = remove_overlapping_rows(unique_transposon_df)
        else:
            non_overlapping_transposon_df = pd.DataFrame()
            print("No matching transposon entries found.")
    else:
        non_overlapping_transposon_df = pd.DataFrame()
    non_overlapping_transposon_df['Category']='Mobile Genetic Element'
    non_overlapping_transposon_df= non_overlapping_transposon_df.rename(columns={"Query Name": "Name of Query"})
    return non_overlapping_transposon_df



def blast_against_replicon_database(replicon_database, plasmid, path):
        #plasmidfinder run 
    if not os.path.exists("makedb_folder/repdb"):
        os.makedirs("makedb_folder/repdb")
        # first make oric database from the fasta file we have. 
    os.system(f'makeblastdb -in {replicon_database} -dbtype nucl -out makedb_folder/repdb/plasmidfinder_database > /dev/null 2>&1')
    # run blast against that oric databse
    os.system(f'blastn -query {path}/{plasmid}.fasta -db makedb_folder/repdb/plasmidfinder_database  -outfmt 5 -out tmp_files/resultplasmidfinder.xml')

    #make a dataframe from the oric balst result
    import pandas as pd
    from Bio.Blast import NCBIXML

    # Initialize the lists to store information
    name_of_query = []
    titlelist = []
    lengthlist = []
    scorelist = []
    gaplist = []
    evallist = []
    sequencelist = []
    startpos = []
    endpos = []
    product_list = []

    # Open and parse the BLAST XML result file
    with open("tmp_files/resultplasmidfinder.xml", "r") as result:
        records = NCBIXML.parse(result)
        item = next(records)
        j = 1  # Start the counter from 1

        for alignment in item.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01 and hsp.identities / alignment.length * 100 > 90:
                    name_of_query.append("Replicon 1" + str(j))
                    startpos.append(hsp.query_start)
                    endpos.append(hsp.query_end)
                    titlelist.append(alignment.title.split(' ')[1])
                    lengthlist.append(alignment.length)
                    scorelist.append(hsp.score)
                    gaplist.append(hsp.gaps)
                    evallist.append(hsp.expect)
                    sequencelist.append(hsp.query)
                    product_list.append('Replicon type')
                    j += 1  # Increment the counter after appending

    # Create the DataFrame from the collected data
    replicon_dataframe = pd.DataFrame({
        "Query Name": name_of_query,
        "Gene Name": titlelist,
        "Length": lengthlist,
        "Score": scorelist,
        "Gaps": gaplist,
        "E-value": evallist,
        "Sequence": sequencelist,
        "Start Position": startpos,
        "End Position": endpos,
        "Product": product_list
    })

    def remove_similar_positions(df):
        df['Start Position Range'] = df['Start Position'] // 50
        df['End Position Range'] = df['End Position'] // 50
        
        # Group by Start Position Range and End Position Range and select row with minimum E-value in each group
        df = df.loc[df.groupby(['Start Position Range', 'End Position Range'])['E-value'].idxmin()]

        # Drop the temporary columns used for grouping
        df = df.drop(columns=['Start Position Range', 'End Position Range'])
        
        return df

    replicon_dataframe = remove_similar_positions(replicon_dataframe)
    replicon_dataframe .reset_index(drop=True, inplace=True)
    replicon_dataframe = replicon_dataframe .rename(columns={"Query Name": "Name of Query"})
    replicon_dataframe ['Category']='Replicon'
    return replicon_dataframe




def process_final_dataframe(selected_rows_per_query, replicon_dataframe, oric_dataframe, orit_dataframe, transposon):
    merged_df = pd.concat([selected_rows_per_query, replicon_dataframe, oric_dataframe, orit_dataframe, transposon], ignore_index=True)
    columns = ['index']
    merged_df.drop(columns, inplace=True, axis=1)
    df_dedup = merged_df.loc[merged_df.groupby('Name of Query')['E-value'].idxmin()]

    fasta_data = {
        "Name of Query": [],
        "Sequence": []
    }

    # Open and read the FASTA file
    with open('tmp_files/Query_Fasta.fsa', 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line.strip().lstrip('>')
                fasta_data["Name of Query"].append(header)
                fasta_data["Sequence"].append(next(file).strip())

    # Convert the dictionary to a DataFrame
    fasta_df = pd.DataFrame(fasta_data)

    # Identify headers from FASTA not present in the original DataFrame
    missing_headers_df = fasta_df[~fasta_df['Name of Query'].isin(df_dedup['Name of Query'])].copy()

    # Read the positions from the text file
    with open('tmp_files/positions.txt', 'r') as file:
        positions = [int(line.strip()) for line in file.readlines()]

    # Extract starting and ending positions
    start_positions = positions[0::2]  # Odd rows: starting positions
    end_positions = positions[1::2]    # Even rows: ending positions

    # Assign start and end positions to the missing headers DataFrame using loc
    indices = missing_headers_df.index
    missing_headers_df.loc[indices, 'Start Position'] = start_positions[:len(missing_headers_df)]
    missing_headers_df.loc[indices, 'End Position'] = end_positions[:len(missing_headers_df)]
    missing_headers_df.loc[indices, 'Length'] = 0
    missing_headers_df.loc[indices, 'Score'] = 0
    missing_headers_df.loc[indices, 'Gaps'] = 0
    missing_headers_df.loc[indices, 'E-value'] = 0
    missing_headers_df.loc[indices, 'Gene Name'] = 'ORF'
    missing_headers_df.loc[indices, 'Category'] = 'Open Reading Frame'
    missing_headers_df.loc[indices, 'Product'] = 'Open Reading Frame'

    # Concatenate the new rows to the original deduplicated DataFrame
    updated_df = pd.concat([df_dedup, missing_headers_df], ignore_index=True)

    # Ensure 'Start Position' is integer
    updated_df['Start Position'] = updated_df['Start Position'].astype(int)

    # Sort the DataFrame by 'Start Position'
    sorted_df = updated_df.sort_values(by='Start Position')
    sorted_df = sorted_df.sort_values(by=['Start Position', 'End Position'])
    sorted_df['Start Position'] = sorted_df['Start Position'].astype(int)
    sorted_df['End Position'] = sorted_df['End Position'].astype(int)

    # Identify rows with overlapping positions
    overlaps = sorted_df.duplicated(subset=['Start Position', 'End Position'], keep=False)

    # Remove rows where Gene Name is "ORF" and there is an overlap
    final_df = sorted_df[~((sorted_df['Gene Name'] == 'ORF') & overlaps)]
    return final_df


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd

def make_genbank_file(only_sequence, final_df, genbank_path, plasmid):
    # Read complement start and end positions from the text files
    with open('tmp_files/complementstarting.txt', 'r') as file:
        complement_starts = [int(line.strip()) for line in file.readlines()]

    with open('tmp_files/complementending.txt', 'r') as file:
        complement_ends = [int(line.strip()) for line in file.readlines()]

    # Handle None sequence by providing an empty string
    sequence_data = '' if only_sequence is None else only_sequence

    record = SeqRecord(
        Seq(sequence_data),  # Initialize with a sequence or an empty string
        id="00000",
        name=plasmid,
        description="Genbank file for Plasmid " + plasmid,
        annotations={"molecule_type": "DNA"}  # Specify the molecule type
    )

    final_df['Start Position'] = pd.to_numeric(final_df['Start Position'], errors='coerce').fillna(0).astype(int)
    final_df['End Position'] = pd.to_numeric(final_df['End Position'], errors='coerce').fillna(0).astype(int)

    # Function to determine if a position is complement
    def is_complement(start, end, comp_starts, comp_ends):
        return start in comp_starts and end in comp_ends

    def determine_feature_type(category):
        if category == "Mobile Genetic Element":
            return "MGE"
        elif category == "Origin of Transfer":
            return "OriT"
        elif category == "Origin of Replication":
            return "OriC"
        elif category == "Replicon":
            return "Replicon"
        else:
            return "CDS"  # Default feature type

    # Add features to the GenBank record
    for idx, row in final_df.iterrows():
        feature_type = determine_feature_type(row['Category'])

        start_position = row['Start Position'] - 1
        end_position = row['End Position']

        # Check if the end position is greater than or equal to the start position
        if end_position < start_position:
            print(f"Skipping invalid feature with start {start_position + 1} and end {end_position}")
            continue

        # Determine if the feature is on the complement strand
        if is_complement(row['Start Position'], row['End Position'], complement_starts, complement_ends):
            location = FeatureLocation(start_position, end_position, strand=-1)
        else:
            location = FeatureLocation(start_position, end_position, strand=1)

        feature = SeqFeature(
            location=location,
            type=feature_type,
            qualifiers={
                "gene": row['Gene Name'],
                "product": row['Product'],
                "translation": row['Sequence'],
                "category": row['Category']
            }
        )
        record.features.append(feature)

    # Save the GenBank file
    with open(genbank_path, 'w') as output_file:
        SeqIO.write(record, output_file, "genbank")

    print(f"GenBank file has been saved to {genbank_path}")


'''def make_genbank_file(only_sequence,final_df,genbank_path,plasmid):
    # Read complement start and end positions from the text files
    with open('tmp_files/complementstarting.txt', 'r') as file:
        complement_starts = [int(line.strip()) for line in file.readlines()]

    with open('tmp_files/complementending.txt', 'r') as file:
        complement_ends = [int(line.strip()) for line in file.readlines()]
    record = SeqRecord(
        Seq(only_sequence),  # Initialize with an empty sequence
        id="00000",
        name=plasmid,
        description="Genbank file for Plasmid"+plasmid,
        annotations={"molecule_type": "DNA"}  # Specify the molecule type
    )

    final_df['Start Position'] = pd.to_numeric(final_df['Start Position'], errors='coerce').fillna(0).astype(int)
    final_df['End Position'] = pd.to_numeric(final_df['End Position'], errors='coerce').fillna(0).astype(int)

    # Function to determine if a position is complement
    def is_complement(start, end, comp_starts, comp_ends):
        return start in comp_starts and end in comp_ends

    def determine_feature_type(category):
        if category == "Mobile Genetic Element":
            return "MGE"
        elif category == "Origin of Transfer":
            return "OriT"
        elif category == "Origin of Replication":
            return "OriC"
        elif category == "Replicon":
            return "Replicon"
        else:
            return "CDS"  # Default feature type

    # Add features to the GenBank record
    for idx, row in final_df.iterrows():
        feature_type = determine_feature_type(row['Category'])
        
        # Determine if the feature is on the complement strand
        if is_complement(row['Start Position'], row['End Position'], complement_starts, complement_ends):
            location = FeatureLocation(row['Start Position']-1, row['End Position'], strand=-1)
        else:
            location = FeatureLocation(row['Start Position']-1, row['End Position'], strand=1)
        
        feature = SeqFeature(
            location=location,
            type=feature_type,
            qualifiers={
                "gene": row['Gene Name'],
                "product": row['Product'],
                "translation":row['Sequence'],
                "category":row['Category']
            }
        )
        record.features.append(feature)

    # Save the GenBank file
    with open(genbank_path, 'w') as output_file:
        SeqIO.write(record, output_file, "genbank")

    print(f"GenBank file has been saved to {genbank_path}")'''



def extract_genbank_info(plasmid,pathofdir):
    file_path = pathofdir+"/"+plasmid+'.gbk'
    record = SeqIO.read(file_path, "genbank")
    
    # Extract CDS translations
    cds_translations = []
    for feature in record.features:
        if feature.type == "CDS":
            translation = feature.qualifiers.get("translation", [""])[0]
            cds_translations.append(translation)
    
    # Extract the whole nucleotide sequence and its length
    whole_nucleotide_sequence = str(record.seq)
    whole_sequence_length = len(record.seq)
    
    return cds_translations, whole_nucleotide_sequence,whole_sequence_length


from Bio import SeqIO, Entrez
from Bio.Seq import UndefinedSequenceError

Entrez.email = "hislam2@ur.rochester.edu"  # Replace with your actual email

def fetch_sequence_from_accession(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq), len(record.seq)
    except Exception as e:
        print(f"Error fetching sequence for {accession_id}: {e}")
        return "", 0

def get_length_from_metadata(record):
    try:
        contig_info = record.annotations["contig"]
        return int(contig_info.split("..")[-1])
    except (KeyError, ValueError):
        return "Unknown"

def extract_genbank_info(plasmid, pathofdir):
    file_path = pathofdir + "/" + plasmid + '.gbk'
    record = SeqIO.read(file_path, "genbank")
    
    # Extract CDS translations
    cds_translations = []
    for feature in record.features:
        if feature.type == "CDS":
            translation = feature.qualifiers.get("translation", [""])[0]
            cds_translations.append(translation)
    
    # Extract the whole nucleotide sequence and its length
    try:
        whole_nucleotide_sequence = str(record.seq)
        whole_sequence_length = len(record.seq)
    except UndefinedSequenceError:
        print(f"Sequence is undefined in the GenBank file: {file_path}. Fetching from NCBI using accession ID.")
        accession_id = record.id
        whole_nucleotide_sequence, whole_sequence_length = fetch_sequence_from_accession(accession_id)
        
        if not whole_nucleotide_sequence:
            whole_sequence_length = get_length_from_metadata(record)
            whole_nucleotide_sequence = ""

    return cds_translations, whole_nucleotide_sequence, whole_sequence_length


def editing_cds_list(CDS_list):
    filtered_list = [item for item in CDS_list if item]
    return filtered_list

import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

def clean_sequence(sequence):
    # Convert to uppercase and remove any character that is not a valid DNA base
    return re.sub(r'[^ATCGN]', '', sequence.upper())

def save_fasta_file(plasmid, nucleotide_sequence):
    # Clean the nucleotide sequence
    cleaned_sequence = clean_sequence(nucleotide_sequence)
    
    # Create the SeqRecord object
    record = SeqRecord(Seq(cleaned_sequence), id=plasmid, description="")
    
    # Define the output file path
    fasta_file_path = os.path.join('tmp_files', plasmid + '.fasta')
    
    try:
        # Write the SeqRecord to a FASTA file
        with open(fasta_file_path, 'w') as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        print(f"FASTA file successfully created: {fasta_file_path}")
    except Exception as e:
        print(f"Failed to write FASTA file: {e}")


def getpositionsofCDS_genbank(plasmid, path):
    file_path = os.path.join(path, f"{plasmid}.gbk")
    list_of_positions = []
    
    try:
        record = SeqIO.read(file_path, "genbank")
        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start) + 1  # Convert to 1-based indexing
                end = int(feature.location.end)
                if start > 0:  # Ensure start is greater than 0
                    list_of_positions.append((start, end))
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return list_of_positions


def complementpositions_genbank(plasmid, path):
    file_path = os.path.join(path, f"{plasmid}.gbk")
    try:
        record = SeqIO.read(file_path, "genbank")
    except FileNotFoundError:
        print(f"Error: The file {file_path} does not exist.")
        return [], []
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return [], []

    list_of_complementstart = []
    list_of_complementend = []

    for feature in record.features:
        if feature.type == "CDS" and feature.location.strand == -1:  # Check if feature is CDS and on the complement strand
            start = int(feature.location.start) + 1  # Adjust to 1-based indexing if needed
            end = int(feature.location.end)
            list_of_complementstart.append(str(start))
            list_of_complementend.append(str(end))
    
    # Ensure the directory exists before trying to write files
    os.makedirs("tmp_files", exist_ok=True)

    # Save the complement start positions to a text file
    with open("tmp_files/complementstarting.txt", "w") as start_file:
        start_file.write("\n".join(list_of_complementstart))

    # Save the complement end positions to a text file
    with open("tmp_files/complementending.txt", "w") as end_file:
        end_file.write("\n".join(list_of_complementend))
    
    return list_of_complementstart, list_of_complementend


    


import os
from Bio.Blast import NCBIXML
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

def run_blast_genbank(sequence, index, list_of_positions):
    query_path = f'tmp_files/query_{index}.fsa'
    result_path = f'tmp_files/result_{index}.xml'
    
    # Write the sequence to a temporary FASTA file
    with open(query_path, 'w') as query_file:
        query_file.write(f'>Coding Sequence{index}\n{sequence}')

    # Check if the sequence is actually written to the file
    if os.stat(query_path).st_size == 0:
        #print(f'Query file {query_path} is empty. No sequence data available for BLAST.')
        return []
    
    # Run BLAST
    command = f'blastp -query {query_path} -db makedb_folder/blastdb/custom_database -outfmt 5 -out {result_path}'
    exit_code = os.system(command)
    if exit_code != 0:
        #print(f'BLAST command failed for query file {query_path}. Exit code: {exit_code}')
        return []

    # Parse BLAST result
    try:
        with open(result_path, 'r') as result_file:
            records = NCBIXML.parse(result_file)
            item = next(records)
            results = []
            for alignment in item.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.001 and hsp.identities / alignment.length * 100 > 60:
                        results.append({
                            "Name of Query": f"Coding Sequence{index}",
                            "Start Position": list_of_positions[index][0],
                            "End Position": list_of_positions[index][1],
                            "Title": alignment.title.split(' ')[1],
                            "Length": alignment.length,
                            "Score": hsp.score,
                            "Pident": hsp.identities / alignment.length * 100,
                            "Gaps": hsp.gaps,
                            "E-value": hsp.expect,
                            "Sequence": hsp.query,
                            "Query Length": len(sequence),
                            "Subject Length": alignment.length
                        })
            #if not results:
                #print(f'No valid BLAST hits for query file {query_path}.')
            return results
    except Exception as e:
        #print(f'Failed to parse BLAST results for {query_path}: {e}')
        return []

def select_gene_row_genbank(df):
    # First try to filter rows with 'Percent Identity' = 100
    filtered_df = df[df['Pident'] >= 75]
    
    # If no rows with 100% identity, or less than 10 rows with 100% identity, consider rows with 'Percent Identity' > 90%
    #if filtered_df.empty or len(filtered_df) < 5:
    if filtered_df.empty:
        filtered_df = df[df['Pident'] > 60]
    
    # Count the occurrences of each gene
    gene_counts = filtered_df['Gene Name'].value_counts()
    max_occurrence = gene_counts.max()
    # Identify genes with the maximum occurrence
    most_frequent_genes = gene_counts[gene_counts == max_occurrence].index.tolist()
    # Filter for the most frequent genes
    most_frequent_df = filtered_df[filtered_df['Gene Name'].isin(most_frequent_genes)]
    # Select the row with the lowest E-value for each gene
    selected_rows = most_frequent_df.loc[most_frequent_df.groupby('Gene Name')['E-value'].idxmin()]
    return selected_rows



def initial_blast_against_database_genbank(list_of_cds, list_of_positions, database_path):
    # Assuming database is loaded from a CSV or similar file
    database = pd.read_csv(database_path)

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_blast_genbank, seq, idx, list_of_positions) for idx, seq in enumerate(list_of_cds)]
        results = [future.result() for future in futures]

    # Flatten list of lists
    flat_results = [item for sublist in results for item in sublist]
    
    # Create DataFrame from results
    final_dataframe = pd.DataFrame(flat_results)
    if not final_dataframe.empty:
        final_dataframe['Gene Name'], final_dataframe['index'] = zip(*final_dataframe['Title'].apply(extract_gene_info))
        final_dataframe.dropna(subset=['Gene Name', 'index'], inplace=True)
        
        # Further processing
        unique_query_names = final_dataframe['Name of Query'].unique()
        selected_rows_per_query = pd.concat([select_gene_row_genbank(final_dataframe[final_dataframe['Name of Query'] == query_name]) for query_name in unique_query_names])
        selected_rows_per_query.reset_index(drop=True, inplace=True)
        
        # Append product and category information
        # Ensure index is integer and within the range
        selected_rows_per_query['Category'] = [database['Category'].iloc[int(index)] if int(index) < len(database) else None for index in selected_rows_per_query['index']]
        selected_rows_per_query['Product'] = [database['Product'].iloc[int(index)] if int(index) < len(database) else None for index in selected_rows_per_query['index']]
        
        # Clean up DataFrame
        selected_rows_per_query.drop(['index', 'Title', 'Pident', 'Query Length', 'Subject Length'], axis=1, inplace=True)

    return selected_rows_per_query

def process_final_dataframe_genbank(selected_rows_per_query, replicon_dataframe, oric_dataframe, orit_dataframe, transposon):
    merged_df = pd.concat([selected_rows_per_query, replicon_dataframe, oric_dataframe, orit_dataframe, transposon], ignore_index=True)
    columns = ['index']
    merged_df.drop(columns, inplace=True, axis=1)
    df_dedup = merged_df.loc[merged_df.groupby('Name of Query')['E-value'].idxmin()]
    df_dedup['Start Position'] = df_dedup['Start Position'].astype(int)
    df_dedup['Endt Position'] = df_dedup['End Position'].astype(int)
    sorted_df = df_dedup.sort_values(by='Start Position')
    return sorted_df

