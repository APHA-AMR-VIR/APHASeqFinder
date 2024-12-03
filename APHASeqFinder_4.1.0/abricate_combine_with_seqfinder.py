#!/usr/bin/env python2


'''
APHASeqfinder
version 4.1.0
submitted to github on 23/12/2021
Javier Nunez, AMR Team, Bacteriology (originally from Nicholas Duggett)
Animal and Plant Health Agency
This script takes abricate output, filters it, takes output, compares it seqfinder output and puts into a new spreadsheet
'''


# Authors : Nicholas Duggett <Nick.Duggett@apha.gov.uk>
# Created : 01/04/19
import sys
import pandas as pd
import numpy as np          
import warnings
import re

# Ignore the specific UserWarning
warnings.filterwarnings("ignore", message="DataFrame columns are not unique")
warnings.filterwarnings("ignore", message="Boolean Series key will be reindexed to match DataFrame index")
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)
#pd.set_option('display.max_colwidth', None)
pd.options.mode.chained_assignment = None

version="4.1.0"

### loading seqfinder and abricate results tables
if len(sys.argv)>1:
    df_abr = pd.read_csv(sys.argv[1], sep='\t')
    df_sf = pd.read_csv(sys.argv[2]) 
    df_original = pd.read_csv(sys.argv[3])
    efsa_dict = pd.read_csv(sys.argv[4])
    first_line = sys.argv[5]
    abricate_version = sys.argv[6]
    mlst_fasta = sys.argv[7]
else:  # just for developing code
    pd.set_option('display.max_columns', None)
    #df_abr = pd.read_csv('/home/nickduggett/mnt/fsx-033/javi_temp/VM0533F/April_2024/seqfinder/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20240322/LREC8042/LREC8042.abricate', sep='\t')
    #df_sf = pd.read_csv('/home/nickduggett/mnt/fsx-033/javi_temp/VM0533F/April_2024/seqfinder/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20240322/LREC8042/LREC8042_CompareTo_AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6_good_snps.csv')
    #df_original= pd.read_csv('/home/nickduggett/mnt/fsx-033/javi_temp/VM0533F/April_2024/seqfinder/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20240322/LREC8042/LREC8042_CompareTo_AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6.csv')
    df_abr = pd.read_csv('/home/nickduggett/mnt/fsx-047/ED1200/Klebsiella/Analysis/seqfinder_outputs/Seqfinder_BL959-1007_Mar2024/Klebsiella_vir_cap_genes_20210105b/20240520/BL1003/BL1003.abricate', sep='\t')
    df_sf = pd.read_csv('/home/nickduggett/mnt/fsx-047/ED1200/Klebsiella/Analysis/seqfinder_outputs/Seqfinder_BL959-1007_Mar2024/Klebsiella_vir_cap_genes_20210105b/20240520/BL1003/BL1003_CompareTo_Klebsiella_vir_cap_genes_20210105b_good_snps.csv')
    df_original= pd.read_csv('/home/nickduggett/mnt/fsx-047/ED1200/Klebsiella/Analysis/seqfinder_outputs/Seqfinder_BL959-1007_Mar2024/Klebsiella_vir_cap_genes_20210105b/20240520/BL1003/BL1003_CompareTo_Klebsiella_vir_cap_genes_20210105b.csv')
    efsa_dict=pd.read_csv('/home/nickduggett/mnt/fsx-057/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.1.0/EFSA_panel/EFSA_antimcriobial_panel_dictionary_20230526.csv')
    abricate_version="1.0.1"
    mlst_fasta= '/home/nickduggett/mnt/fsx-057/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.1.0/references/ECO-MLST-MG1655-alleles.fna'
    with open('/home/nickduggett/seqfinder_testing/PSBM-EX017.fasta') as f:
        first_line = f.readline().strip()

df_sf = df_sf.drop(['warnings'], axis=1)
      
# Determining the versions of the assemblers used
match = re.search(r'Unicycler.*', first_line)
if match:
    assembler_versions = match.group()
else:
    assembler_versions = "Unknown"

# Extract the unicycler version number
match = re.search(r'v([\d.]+)', assembler_versions)
unicycler_version = match.group(1) if match else "Unknown"

# Extract the spades version number
match = re.search(r'SPAdes_v([\d.]+)', assembler_versions)
spades_version = match.group(1) if match else "Unknown"

###grouping plasmids by contig and choosing a representative 
#Make list of genes from the SeqFinder output
df_abr["LOCATION"]=df_abr['START'].astype(str)+".."+df_abr['END'].astype(str)
'''
plasmids_only=df_abr[(df_abr.GENE.str.contains("PlasR"))]


# Create a DataFrame to store the dropped rows
dropped_rows = df_abr[~df_abr.index.isin(plasmids_only.index)]

grouped_plasmids = plasmids_only.loc[plasmids_only.groupby("SEQUENCE")["%COVERAGE"].idxmax()]

# Identify the rows that were not included in the grouped part
not_grouped_rows = plasmids_only[~plasmids_only.index.isin(grouped_plasmids.index)]

# Identify the rows that were not included in the grouped part
# Display the plasmids_only DataFrame and the grouped_plasmids DataFrame
print("\nNot Grouped Rows:")
print(not_grouped_rows)



#### removed indexes with plasmids and merging the rest with the representatives
drop_plasmids=df_abr[(~df_abr.GENE.str.contains("PlasR"))]
merge_grouped_plasmids=drop_plasmids,grouped_plasmids
df_abr=pd.concat(merge_grouped_plasmids)
'''

def extract_coverage_info(coverage_str):
    coverage_parts = coverage_str.split("/")
    coverage_start, coverage_end = map(int, coverage_parts[0].split("-"))
    gene_length = int(coverage_parts[1])
    length_of_hit=(coverage_end-coverage_start+1)
    return coverage_start, coverage_end, gene_length, length_of_hit

df_abr["COVERAGE_START"], df_abr["COVERAGE_END"], df_abr["GENE_LENGTH"],df_abr["LENGTH_OF_HIT"] = zip(*df_abr["COVERAGE"].apply(extract_coverage_info))

# Calculate the calculated identity based on the proportion identity and %IDENTITY
df_abr["CALCULATED_IDENTITY"] = df_abr["LENGTH_OF_HIT"] * (df_abr["%IDENTITY"] / 100)

def combine_rows(group):
    if group["GENE"].nunique() == 1:  # Check if there are multiple identical strings
        min_coverage_start = group["COVERAGE_START"].min()
        max_coverage_end = group["COVERAGE_END"].max()
        gene_length = group["GENE_LENGTH"].iloc[0]

        required_coverage = gene_length * 0.9  # Calculate the required 90% coverage

        if max_coverage_end - min_coverage_start + 1 >= required_coverage and group["%COVERAGE"].max() < 90:
            # If the combined range covers at least 90% of the gene length, return the combined row
            combined_row = group.iloc[[0]].copy()
            combined_row["SEQUENCE"] = ";".join(group["SEQUENCE"])
            combined_coverage = ";".join(f"{start}-{end}" for start, end in zip(group["COVERAGE_START"], group["COVERAGE_END"]))
            combined_row["COVERAGE"] = f"{combined_coverage}/{gene_length}"
            combined_row["%COVERAGE"] = ((max_coverage_end - min_coverage_start + 1) / gene_length) * 100

            # Concatenate the "LOCATION" column values with a semicolon
            combined_row["LOCATION"] = ";".join(group["LOCATION"])

            # Calculate the sum of "CALCULATED_IDENTITY" and "LENGTH_OF_HIT" for the combined rows
            sum_calculated_identity = group["CALCULATED_IDENTITY"].sum()
            sum_length_of_hit = group["LENGTH_OF_HIT"].sum()

            # Calculate the new "%IDENTITY" for the combined rows
            new_identity = sum_calculated_identity / sum_length_of_hit
            combined_row["%IDENTITY"] = round(new_identity * 100, 2)

            return combined_row
        else:
            # If the combined range does not cover at least 90% of the gene length, return the original rows
            return group.copy()
    else:
        # If there are no multiple identical strings, return the original rows
        return group.copy()
    
# Apply the combine_rows function to each group (grouped by "X.FILE" and "GENE")
df_abr = df_abr.groupby(["#FILE", "GENE"], group_keys=True).apply(combine_rows).reset_index(drop=True)

def extract_string_after_underscore_abr(row):
    if "PlasR" in row["GENE"]:
        plasmid_id = row["GENE"].split("_", 1)[-1]
        return plasmid_id
    return ""

df_abr["plasmid_result_abr"] = df_abr.apply(extract_string_after_underscore_abr, axis=1)


#Circumvents the filtering rule that would exclude ant3, gyrA, parC and or ampP
circumvent_rule = df_abr[((df_abr.GENE.str.contains("ant3")) & (df_abr['%COVERAGE'] >= 80)) | (df_abr.GENE.str.contains("gyrA")) | (df_abr.GENE.str.contains("parC")) | (df_abr.GENE.str.contains("ampP"))]
#Drop hits for chromosomal genes from the abricate output, ones less than 90% coverage and identity, and those that are plasmid genes
filtering_rule = df_abr[(~df_abr.GENE.str.contains("Chr") &(df_abr['%COVERAGE'] >= 80) &(df_abr['%IDENTITY'] >= 80) )]
# Exception for the rows containing "gyrA", "parC" or "ampP"
exception_conditions = ((df_abr.GENE.isin(["gyrA", "parC", "ampP"])))
#Merge ant3 back into the dataframe after circumventing the filter
merge_ant3 = filtering_rule,circumvent_rule
filtering_rule = pd.concat(merge_ant3)
#Use same rule as SeqFinder to make a 4 character gene code from the abricate hits and place into a new column called GENE_id
df_abr_GENE = filtering_rule['GENE']

gene_prefixes=list()
for i in df_abr_GENE:
    s = i.find('_')
    gene_prefixes.append(i[s+1:s+5])
filtering_rule['GENE_id']=gene_prefixes


#Filter the GENE_id list with the ones found in the SeqFinder dataframe
df_sf_genes= (df_sf['gene'])
filter_for_sf_genes = filtering_rule.loc[filtering_rule['GENE_id'].isin(df_sf_genes)]
GENE_id = (filter_for_sf_genes['GENE_id'])

#Create a dictionary of filtered GENE_id hits and the GENE hit
grouped_df_gene= filter_for_sf_genes.groupby('GENE_id')['GENE'].apply(','.join).reset_index()
abricate_dictionary_gene = grouped_df_gene.set_index('GENE_id').T.to_dict('list')

abricate_dictionary_sequence = filter_for_sf_genes[['GENE','SEQUENCE']]
abricate_dictionary_sequence = abricate_dictionary_sequence.set_index('GENE').T.to_dict('list')

abricate_dictionary_coverage = filter_for_sf_genes[['GENE','%COVERAGE']]
abricate_dictionary_coverage = abricate_dictionary_coverage.set_index('GENE').T.to_dict('list')

abricate_dictionary_coverage_map = filter_for_sf_genes[['GENE','COVERAGE']]
abricate_dictionary_coverage_map = abricate_dictionary_coverage_map.set_index('GENE').T.to_dict('list')

abricate_dictionary_identity = filter_for_sf_genes[['GENE','%IDENTITY']]
abricate_dictionary_identity = abricate_dictionary_identity.set_index('GENE').T.to_dict('list')

abricate_dictionary_location = filter_for_sf_genes[['GENE','LOCATION']]
abricate_dictionary_location = abricate_dictionary_location.set_index('GENE').T.to_dict('list')


def clean_df_sf(df_sf, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_identity, abricate_dictionary_location, abricate_dictionary_coverage, abricate_dictionary_coverage_map):
    # Map the dictionaries onto the SeqFinder dataframe into columns listed below
    df_sf['result_abr'] = df_sf['gene'].map(abricate_dictionary_gene)
    df_sf['contig_abr'] = df_sf['id'].map(abricate_dictionary_sequence)
    df_sf['identity_abr'] = df_sf['id'].map(abricate_dictionary_identity)
    df_sf['location_abr'] = df_sf['id'].map(abricate_dictionary_location)
    df_sf['coverage_abr'] = df_sf['id'].map(abricate_dictionary_coverage)
    df_sf['coverage_map_abr'] = df_sf['id'].map(abricate_dictionary_coverage_map)

    # Clean up the abricate columns to only include the text and no brackets
    df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[2:-2])
    df_sf['contig_abr'] = df_sf['contig_abr'].map(lambda x: str(x)[2:-2])
    df_sf['identity_abr'] = df_sf['identity_abr'].map(lambda x: str(x)[1:-1])
    df_sf['location_abr'] = df_sf['location_abr'].map(lambda x: str(x)[2:-2])
    df_sf['coverage_abr'] = df_sf['coverage_abr'].map(lambda x: str(x)[1:-1])
    df_sf['coverage_map_abr'] = df_sf['coverage_map_abr'].map(lambda x: str(x)[2:-2])

    # Replace 'a' (a carry over from NaN) with blanks
    df_sf = df_sf.replace({'a': ''})
    
    return df_sf
df_sf = clean_df_sf(df_sf, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_identity, abricate_dictionary_location, abricate_dictionary_coverage, abricate_dictionary_coverage_map)

#print(df_sf)
if "class" in df_sf.columns:
    df_sf_plasmids = df_sf[(df_sf['class'].str.contains("PlasR")) & (df_sf['coverage_abr'] != "")]
#print(df_sf_plasmids)

def extract_string_after_underscore(row):
    if "PlasR" in row["class"]:
        plasmid_id = row["id"].split("_", 1)[-1]
        plasmid_result_abr = row["result_abr"].split("_", 1)[-1]
        return plasmid_id, plasmid_result_abr
    return "", ""

if 'EFSA_dict' in df_sf.columns:
    df_sf[["plasmid_result_seq", "plasmid_result_abr"]] = df_sf.apply(extract_string_after_underscore, axis=1, result_type="expand")

    condition = (df_sf["plasmid_result_seq"] == df_sf["plasmid_result_abr"]) & (df_sf["plasmid_result_seq"] != "") & (df_sf["contig_abr"] == "")
    
    if condition.any():
        # Convert the "contig_abr" column to strings if it contains lists
        df_sf.loc[condition, 'contig_abr'] = df_sf.loc[condition, 'contig_abr'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)
    
        abricate_dictionary_sequence_plasmid = filter_for_sf_genes[['plasmid_result_abr','SEQUENCE']]
        abricate_dictionary_sequence_plasmid = abricate_dictionary_sequence_plasmid.set_index('plasmid_result_abr').T.to_dict('list')
    
        abricate_dictionary_identity_plasmid = filter_for_sf_genes[['plasmid_result_abr','%IDENTITY']]
        abricate_dictionary_identity_plasmid = abricate_dictionary_identity_plasmid.set_index('plasmid_result_abr').T.to_dict('list')
    
        abricate_dictionary_location_plasmid = filter_for_sf_genes[['plasmid_result_abr','LOCATION']]
        abricate_dictionary_location_plasmid = abricate_dictionary_location_plasmid.set_index('plasmid_result_abr').T.to_dict('list')
        
        abricate_dictionary_coverage_plasmid = filter_for_sf_genes[['plasmid_result_abr','%COVERAGE']]
        abricate_dictionary_coverage_plasmid = abricate_dictionary_coverage_plasmid.set_index('plasmid_result_abr').T.to_dict('list')
    
        abricate_dictionary_coverage_map_plasmid = filter_for_sf_genes[['plasmid_result_abr','COVERAGE']]
        abricate_dictionary_coverage_map_plasmid = abricate_dictionary_coverage_map_plasmid.set_index('plasmid_result_abr').T.to_dict('list')
    
        df_sf.loc[condition, 'contig_abr'] = df_sf.loc[condition, 'plasmid_result_seq'].map(abricate_dictionary_sequence_plasmid)
        df_sf.loc[condition, 'identity_abr'] = df_sf.loc[condition, 'plasmid_result_seq'].map(abricate_dictionary_identity_plasmid)
        df_sf.loc[condition, 'location_abr'] = df_sf.loc[condition, 'plasmid_result_seq'].map(abricate_dictionary_location_plasmid)
        df_sf.loc[condition, 'coverage_abr'] = df_sf.loc[condition, 'plasmid_result_seq'].map(abricate_dictionary_coverage_plasmid)
        df_sf.loc[condition, 'coverage_map_abr'] = df_sf.loc[condition, 'plasmid_result_seq'].map(abricate_dictionary_coverage_map_plasmid)
        
        df_sf.loc[condition, 'contig_abr'] = df_sf.loc[condition, 'contig_abr'].apply(lambda x: str(x)[2:-2])
        df_sf.loc[condition, 'identity_abr'] = df_sf.loc[condition, 'identity_abr'].apply(lambda x: str(x)[1:-1])
        df_sf.loc[condition, 'location_abr'] = df_sf.loc[condition, 'location_abr'].apply(lambda x: str(x)[2:-2])
        df_sf.loc[condition, 'coverage_abr'] = df_sf.loc[condition, 'coverage_abr'].apply(lambda x: str(x)[1:-1])
        df_sf.loc[condition, 'coverage_map_abr'] = df_sf.loc[condition, 'coverage_map_abr'].apply(lambda x: str(x)[2:-2])
        

#Pick the plasmid rep gene with the highest coverage as the one to report
list_of_contigs = df_sf['contig_abr']
filter_for_abr_contigs = df_abr.loc[df_abr['SEQUENCE'].isin(list_of_contigs)]
filter_for_abr_contigs = filter_for_abr_contigs[(df_abr.GENE.str.contains("PlasR"))]
#grouped_plasmids = df_abr.groupby(by="SEQUENCE")['%COVERAGE'].max()
abricate_dictionary_plasmid = filter_for_abr_contigs[['SEQUENCE','GENE']]
abricate_dictionary_plasmid = abricate_dictionary_plasmid.set_index('SEQUENCE').T.to_dict('list')

def process_plasmid_abr(df_sf, abricate_dictionary_plasmid):
    df_sf['plasmid_abr'] = df_sf['contig_abr'].map(abricate_dictionary_plasmid)
    df_sf['plasmid_abr'] = df_sf['plasmid_abr'].map(lambda x: str(x)[2:-2])
    return df_sf

df_sf = process_plasmid_abr(df_sf, abricate_dictionary_plasmid)

#print(df_sf)
###Temporary workaround to fix the problem of the same gene being present in the sample twice
column_name = "result_abr"
df_sf[column_name] = df_sf[column_name].apply(lambda x: ",".join(pd.Series(x).str.split(",").explode().str.strip().drop_duplicates()))

def process_multiple_strings(row):
    result_abr = row['result_abr']
    id_string = row['id']
    if isinstance(result_abr, str) and ',' in result_abr:
        result_strings = result_abr.split(',')
        matched_string = [string for string in result_strings if string == id_string]
        if matched_string:
            row['result_abr'] = matched_string[0]
        else:
            row['result_abr'] = ''
    return row


if df_sf['result_abr'].str.contains(',').any():
    df_sf['result_abr'] = df_sf['result_abr'].str.split(',')
    df_sf['multiple_strings'] = df_sf['result_abr'].apply(lambda x: len(x) > 1)

    # Create a set of unique ids
    unique_ids = set(df_sf['id'])
    
    column_mapping = {'ref_len':'ref_len', 'mapped_len':'mapped_len', 'mean_depth':'mean_depth',
                      'norm_depth':'norm_depth', 'non_calls':'non_calls','perc_mapped':'perc_mapped',
                      'snps':'snps','other':'other','good_snps':'good_snps','syn':'syn','non':'non',
                      'annotation':'annotation'}
    
    
    
    # Function to create new rows for missing strings
    def create_missing_rows(row):
        missing_strings = [id_str for id_str in row['result_abr'] if id_str not in unique_ids]
        new_rows = []
        for string in missing_strings:
            new_row = {'id': string}
            matching_row = df_sf.loc[(df_sf['result_abr'].apply(lambda x: string in x)) & (df_sf['id'] != string)]
            new_row['strain'] = matching_row['strain'].iloc[0]
            new_row['gene'] = matching_row['gene'].iloc[0]
            if 'class' in df_sf.columns:
                new_row['class'] = matching_row['class'].iloc[0]
            for col, mapping in column_mapping.items():
                new_row[col] = df_original[df_original['id'] == string][mapping].iloc[0]
            new_rows.append(new_row)
            #df_sf.loc[idx, 'result_abr'] = [id_str for id_str in df_sf.loc[idx, 'result_abr'] if id_str != string]
            df_sf.at[idx, 'result_abr'] = [id_str for id_str in df_sf.at[idx, 'result_abr'] if id_str != string]
            new_row['result_abr'] = string
        return new_rows
  
    # Apply the function to rows with multiple strings and missing ids
    new_rows = []
    for idx, row in df_sf.iterrows():
        if row['multiple_strings'] and any(id_str not in unique_ids for id_str in row['result_abr']):
            new_rows.extend(create_missing_rows(row))
    if new_rows:
        # Append new rows to the DataFrame
        #df_sf = df_sf.append(new_rows, ignore_index=True) #Original code
        new_rows_df = pd.DataFrame(new_rows)
    
        # Apply the functions from before that iterate through the dataframe again
        # with the updated information to the new rows that were created to stop the
        # overwriting of plasmids of the same Inc-type being lost
        new_rows_df=clean_df_sf(new_rows_df, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_identity, abricate_dictionary_location, abricate_dictionary_coverage, abricate_dictionary_coverage_map)
        new_rows_df = process_plasmid_abr(new_rows_df, abricate_dictionary_plasmid)
        for index, row in new_rows_df.iterrows():
            if 'PlasR' in row['id']:
                new_rows_df.at[index, 'plasmid_abr'] = row['id']
        df_sf = pd.concat([df_sf, new_rows_df], ignore_index=True)
        # Drop the intermediate column
        df_sf.drop('multiple_strings', axis=1, inplace=True)
        df_sf = df_sf.apply(process_multiple_strings, axis=1)
        # Re-apply the efsa dictionary for the new row
        efsa_as_dict = efsa_dict.set_index('id').T.to_dict('list')
        if "antimicrobial" in df_sf.columns:
            df_sf['antimicrobial']= df_sf['id'].map(efsa_as_dict)
        # Forward fill the empty values in the "reference" column
        df_sf['reference'].fillna(method='ffill', inplace=True)
        # Forward fill the empty values in the "EFSA_dict" column
        if "EFSA_dict" in df_sf.columns:
            df_sf['EFSA_dict'].fillna(method='ffill', inplace=True)
        # Formatting the dataframe so it looks a bit nicer
        df_sf=df_sf.fillna('')
        if "antimicrobial" in df_sf.columns:
            df_sf['antimicrobial'] = df_sf['antimicrobial'].map(lambda x: str(x)[2:-2])
        df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[2:-2])
        #Drop duplicate ids that get included sometimes
        df_sf = df_sf.drop_duplicates(subset='id')
        df_sf = df_sf.replace({'a': ''})
    elif not new_rows:
        df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[1:-1])


# Filtering some pesky genes
values_to_check = ["TEM", "SHV", "CTX", "OXA", "CMY", "oqxA", "oqxB", "tet", "aadA","dfrA","mcr","ant2"]
mask = df_sf["id"].str.contains('|'.join(values_to_check))
blank_result_abr_mask = mask & ~df_sf["result_abr"].str.contains("-")

# Create an additional mask to check if "result_abr" strings are found anywhere in the "id" column
result_abr_found_mask = df_sf["result_abr"].str.contains('|'.join(df_sf["id"].apply(re.escape)))



# Update the original masks based on the "result_abr_found_mask"
mask = mask & result_abr_found_mask
blank_result_abr_mask = blank_result_abr_mask & result_abr_found_mask

if df_sf[mask].shape[0] > 0:
    if mask.any() and ((df_sf.loc[mask, "id"] == df_sf.loc[mask, "result_abr"]).any() or blank_result_abr_mask.any()):
        df_sf = df_sf.loc[(mask & (df_sf["id"] == df_sf["result_abr"])) | ~mask | blank_result_abr_mask]
else:
    df_sf = df_sf
    
#Reformat the output columns in the order we desire
abricate_columns=(df_sf.loc[:,['plasmid_abr','identity_abr','coverage_abr','coverage_map_abr','location_abr']])
df_sf=df_sf.iloc[:,:-5]
merge=df_sf,abricate_columns

df_sf=pd.concat([df_sf, abricate_columns], axis=1, sort=False)

# Adding a new column to the dataframe to tell the user if multiple genes are found on the same contig
df_sf['genes_on_same_contig'] = df_sf[df_sf['contig_abr'].str.strip() != ''].groupby('contig_abr')['id'].transform(lambda x: ', '.join(x) if x.astype(str).str.strip().str.len().notnull().all() and len(x) > 1 else '')

#Put versions of tools used into final output
df_sf['SeqFinder_version']=version

df_sf['Abricate_version']=abricate_version

df_sf['Unicycler_version']=unicycler_version

df_sf['SPAdes_version']=spades_version


#Move columns around to make things look nicer
if 'EFSA_dict' in df_sf.columns:
    cols = df_sf.pop('reference'), df_sf.pop('EFSA_dict')
    df_sf['reference'], df_sf['EFSA_dict'] = cols
elif 'Vir_dict' in df_sf.columns:
    cols = df_sf.pop('reference'), df_sf.pop('Vir_dict')
    df_sf['reference'], df_sf['Vir_dict'] = cols
elif 'Bio_metal_dict' in df_sf.columns:
    cols = df_sf.pop('reference'), df_sf.pop('Bio_metal_dict')
    df_sf['reference'], df_sf['Bio_metal_dict'] = cols 
else:
    cols = df_sf.pop('reference')
    df_sf['reference'] = cols

#Delete uneeded columns from non-AMR dataframes
if 'EFSA_dict' not in df_sf.columns:
    df_sf = df_sf.drop('plasmid_abr', axis=1)

#Put MLST scheme used in final output
df_sf['MLST_reference']=mlst_fasta.split("/")[-1]

#Getting rid of any "NaN"s carried over
df_sf.replace("NaN", np.nan, inplace=True)
df_sf.fillna("", inplace=True)

# Remove duplicate columns
df_sf = df_sf.loc[:, ~df_sf.columns.duplicated()]

# Convert the "coverage_abr" column to numeric (float)
df_sf["coverage_abr"] = pd.to_numeric(df_sf["coverage_abr"], errors="coerce")
# Convert the "identity_abr" column to numeric (float)
df_sf["identity_abr"] = pd.to_numeric(df_sf["identity_abr"], errors="coerce")


def generate_warning(row):
    warnings = []
    if row["perc_mapped"] < 95:
        warnings.append ("Percentage mapped is under 95%")
        
    # Check if the value in the "norm" column is less than 0.4
    if row["norm_depth"] < 0.4:
        warnings.append("Normalised depth is low")

    # Check if there is a ";" present in the "contig_abr" column
    num_contigs = row["contig_abr"].count(";") + 1
    if num_contigs > 1:
        warnings.append(f"Gene split over multiple contigs in abricate ({num_contigs} contigs)")

    # Check if there is a mismatch between "id" and "result_abr" columns
    if row["result_abr"] != "":
        if row["id"] != row["result_abr"]:
            if row["result_abr"] != "" and "PlasR" in row["result_abr"]:
                warnings.append("Same Inc-type but plasmid mismatch between seqfinder and abricate")
            else:
                warnings.append("Mismatch between seqfinder and abricate")
                
    # Check if there is no corresponding hit in abricate
    if row["result_abr"] == "" :
        warnings.append("No corresponding hit in abricate")
    
    if row["coverage_abr"] < 90:
        warnings.append("Abricate gene coverage is below 90%")
        
    if row["identity_abr"] < 90:
        warnings.append("Abricate gene identity is below 90%")
        
    # Concatenate multiple warning messages using a ";"
    return ";".join(warnings)

# Create a new "warnings" column by applying the generate_warning function to each row
df_sf["perc_mapped"] = pd.to_numeric(df_sf["perc_mapped"], errors='coerce')
df_sf["norm_depth"] = pd.to_numeric(df_sf["norm_depth"], errors='coerce')
df_sf["warnings"] = df_sf.apply(generate_warning, axis=1)
# Reorder the columns to make the "warnings" column the first column
df_sf = df_sf[["warnings"] + df_sf.columns[:-1].tolist()]

#Re-introduce any genes that were present twice.
df_sf['result_abr'] = df_sf['gene'].map(abricate_dictionary_gene)
df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[2:-2])
df_sf.replace("NaN", np.nan, inplace=True)
df_sf.fillna("", inplace=True)

#####Make this a mask rather than if rules because it needs to apply per row not per column
# Fixing the ant3-1a/ant3-Ia and aadA1b multiple hit issue
# If the DataFrame contains ant3 and aadA1b, then check further
if df_sf['id'].str.contains('ant3').any() and df_sf['id'].str.contains('aadA1b').any():
    # If the DataFrame contains aadA1b with a warning, then continue
    if df_sf['id'].str.contains('aadA1b').any() and df_sf['warnings'].str.contains("Mismatch between seqfinder and abricate").any():
        # If ant3 is present and the warnings column either contains "Mismatch between seqfinder and abricate" or is blank, then remove the aadA1b row
        if (df_sf['id'].str.contains('ant3')).any() and (((df_sf['warnings'].str.contains("Mismatch between seqfinder and abricate")) | (df_sf['warnings'] == "")).any()):
            df_sf = df_sf[~df_sf['id'].str.contains('aadA1b')]

#####Make this a mask rather than if rules because it needs to apply per row not per column
# Fixing the IncI1 and IncI-G multiple hit issue
# If the DataFrame contains I1 and IncI-G, then check further
if df_sf['id'].str.contains('_I1').any() and df_sf['id'].str.contains('IncI-G').any():
    # If the DataFrame contains IncI-G with a warning about not being found, then continue
    if df_sf['id'].str.contains('IncI-G').any() and df_sf['warnings'].str.contains("No corresponding hit in abricate").any():
        # If I1 is present and the warnings column either contains "Mismatch between seqfinder and abricate" or is blank, then remove the IncI-G row
        if (df_sf['id'].str.contains('I1')).any() and (((df_sf['warnings'].str.contains("Mismatch between seqfinder and abricate")) | (df_sf['warnings'] == "")).any()):
            df_sf = df_sf[~df_sf['id'].str.contains('IncI-G')]
    
    

mask = df_sf['id'].str.contains('ant3') & (df_sf['warnings'].str.contains("Mismatch between seqfinder and abricate")) #& (df_sf['warnings'] != '')
if mask.any():
    df_sf_no_ant3 = df_sf[~df_sf['id'].str.contains('ant3')]
    ant3_df=df_sf[((df_sf.id.str.contains("ant3")))]
    
    #Create dictionaries based on the gene code rather than the full ID as was done previously to ensure any hit for "ant3" is available to be mapped
    abricate_dictionary_sequence = filter_for_sf_genes[['GENE_id','SEQUENCE']]
    abricate_dictionary_sequence = abricate_dictionary_sequence.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_coverage = filter_for_sf_genes[['GENE_id','%COVERAGE']]
    abricate_dictionary_coverage = abricate_dictionary_coverage.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_coverage_map = filter_for_sf_genes[['GENE_id','COVERAGE']]
    abricate_dictionary_coverage_map = abricate_dictionary_coverage_map.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_identity = filter_for_sf_genes[['GENE_id','%IDENTITY']]
    abricate_dictionary_identity = abricate_dictionary_identity.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_location = filter_for_sf_genes[['GENE_id','LOCATION']]
    abricate_dictionary_location = abricate_dictionary_location.set_index('GENE_id').T.to_dict('list')
    
    #Map the dictionaries to the dataframe
    ant3_df['result_abr'] = ant3_df['gene'].map(abricate_dictionary_gene)
    ant3_df['contig_abr'] = ant3_df['gene'].map(abricate_dictionary_sequence)
    ant3_df['identity_abr'] = ant3_df['gene'].map(abricate_dictionary_identity)
    ant3_df['location_abr'] = ant3_df['gene'].map(abricate_dictionary_location)
    ant3_df['coverage_abr'] = ant3_df['gene'].map(abricate_dictionary_coverage)
    ant3_df['coverage_map_abr'] = ant3_df['gene'].map(abricate_dictionary_coverage_map)
    
    #Clean up the re-mapped columns to only include the text and no brackets
    ant3_df['result_abr'] = ant3_df['result_abr'].map(lambda x: str(x)[2:-2])
    ant3_df['contig_abr'] = ant3_df['contig_abr'].map(lambda x: str(x)[2:-2])
    ant3_df['identity_abr'] = ant3_df['identity_abr'].map(lambda x: str(x)[1:-1])
    ant3_df['location_abr'] = ant3_df['location_abr'].map(lambda x: str(x)[2:-2])
    ant3_df['coverage_abr'] = ant3_df['coverage_abr'].map(lambda x: str(x)[1:-1])
    ant3_df['coverage_map_abr'] = ant3_df['coverage_map_abr'].map(lambda x: str(x)[2:-2])

    
    #Add a new warning to the user of the fact ant3 mismatched but was still pulled through
    ant3_df["warnings"] = "Mismatch of ant3 between SeqFinder and Abricate"
    #Merge the dataframes that were split to deal with the ant3 mismatch, also accounting for any "NaN"s that are present
    merge_dfs = df_sf_no_ant3,ant3_df
    df_sf_no_ant3.replace("NaN", np.nan, inplace=True)
    df_sf_no_ant3.fillna("", inplace=True)
    df_sf = pd.concat(merge_dfs)

#Dealing with difference between amino-g0116_ant2-Ia and amino-g2875_ant2-Ia
mask = df_sf['id'].str.contains('ant2') & (df_sf['warnings'].str.contains("Mismatch between seqfinder and abricate")) #& (df_sf['warnings'] != '')
if mask.any():
    df_sf_no_ant2 = df_sf[~df_sf['id'].str.contains('ant2')]
    ant2_df=df_sf[((df_sf.id.str.contains("ant2")))]
    
    #Create dictionaries based on the gene code rather than the full ID as was done previously to ensure any hit for "ant2" is available to be mapped
    abricate_dictionary_sequence = filter_for_sf_genes[['GENE_id','SEQUENCE']]
    abricate_dictionary_sequence = abricate_dictionary_sequence.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_coverage = filter_for_sf_genes[['GENE_id','%COVERAGE']]
    abricate_dictionary_coverage = abricate_dictionary_coverage.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_coverage_map = filter_for_sf_genes[['GENE_id','COVERAGE']]
    abricate_dictionary_coverage_map = abricate_dictionary_coverage_map.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_identity = filter_for_sf_genes[['GENE_id','%IDENTITY']]
    abricate_dictionary_identity = abricate_dictionary_identity.set_index('GENE_id').T.to_dict('list')
    abricate_dictionary_location = filter_for_sf_genes[['GENE_id','LOCATION']]
    abricate_dictionary_location = abricate_dictionary_location.set_index('GENE_id').T.to_dict('list')
    
    #Map the dictionaries to the dataframe
    ant2_df['result_abr'] = ant2_df['gene'].map(abricate_dictionary_gene)
    ant2_df['contig_abr'] = ant2_df['gene'].map(abricate_dictionary_sequence)
    ant2_df['identity_abr'] = ant2_df['gene'].map(abricate_dictionary_identity)
    ant2_df['location_abr'] = ant2_df['gene'].map(abricate_dictionary_location)
    ant2_df['coverage_abr'] = ant2_df['gene'].map(abricate_dictionary_coverage)
    ant2_df['coverage_map_abr'] = ant2_df['gene'].map(abricate_dictionary_coverage_map)
    
    #Clean up the re-mapped columns to only include the text and no brackets
    ant2_df['result_abr'] = ant2_df['result_abr'].map(lambda x: str(x)[2:-2])
    ant2_df['contig_abr'] = ant2_df['contig_abr'].map(lambda x: str(x)[2:-2])
    ant2_df['identity_abr'] = ant2_df['identity_abr'].map(lambda x: str(x)[1:-1])
    ant2_df['location_abr'] = ant2_df['location_abr'].map(lambda x: str(x)[2:-2])
    ant2_df['coverage_abr'] = ant2_df['coverage_abr'].map(lambda x: str(x)[1:-1])
    ant2_df['coverage_map_abr'] = ant2_df['coverage_map_abr'].map(lambda x: str(x)[2:-2])

    
    if ant2_df['id'] != ant2_df['result_abr']:
        #Add a new warning to the user of the fact ant2 mismatched but was still pulled through
        ant2_df["warnings"] = "Mismatch of ant2 between SeqFinder and Abricate"
    else:
        ant2_df["warnings"] =""
    #Merge the dataframes that were split to deal with the ant2 mismatch, also accounting for any "NaN"s that are present
    merge_dfs = df_sf_no_ant2,ant2_df
    df_sf_no_ant2.replace("NaN", np.nan, inplace=True)
    df_sf_no_ant2.fillna("", inplace=True)
    df_sf = pd.concat(merge_dfs)
    
    
#####Make this a mask rather than if rules because it needs to apply per row not per column
#Adding separate warning here for multiple genes being present
'''
if df_sf["result_abr"].str.contains(",").any():
    #If multiple genes detected and the warnings columns is not empty add ";Multiple genes present", otherwise if warnings is empty then just add "Multiple genes present"
    if df_sf['result_abr'].str.contains(',').any() and (df_sf['warnings'] != "").any():
        df_sf.loc[df_sf["result_abr"].str.contains(","), "warnings"] += ";Multiple genes from same family present"
    elif df_sf['result_abr'].str.contains(',').any() and (df_sf['warnings'] == "").any():
        df_sf.loc[df_sf["result_abr"].str.contains(","), "warnings"] += "Multiple genes from same family present"
'''

def update_warnings(row):
    if ',' in row['result_abr']:
        # Split the result_abr values by comma and remove leading/trailing whitespaces
        genes = [gene.strip() for gene in row['result_abr'].split(',')]
        
        # Check for duplicates based on the first 4 characters after an underscore
        seen = set()
        for gene in genes:
            # Extract the substring after the underscore
            substring = gene.split('_')[1][:4]
            
            # Check if the substring is already in the set
            if substring in seen:
                # Update the warnings column if duplicates are found
                if row['warnings'] != "":
                    return row['warnings'] + ";Multiple genes from same family present"
                else:
                    return "Multiple genes from same family present"
            else:
                seen.add(substring)
    
    # If no duplicates are found or if there is no comma in result_abr, return the existing warnings
    return row['warnings']

# Apply the custom function to each row
df_sf['warnings'] = df_sf.apply(update_warnings, axis=1)

#Printing for developing only
if len(sys.argv)>1:
    output_filename = sys.argv[2][:-4]+'_abricate_seqfinder.csv'
    df_sf.to_csv(output_filename,index=False)
    print('Done! Saved output as {}.'.format(output_filename))
else:
    print(df_sf)
    print("DONE!")
   

###To automatically solve the CTX, TEM and CMY multiple occurance issue this can be turned on
#id_list = ['CTX-','TEM-','CMY-']
#filter_for_TEM = df_sf.loc[df_sf['gene'].isin(id_list)]
#grouped_df = filter_for_TEM.groupby(by="gene")['snps'].min()




