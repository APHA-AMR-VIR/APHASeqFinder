#!/usr/bin/env python2


'''
APHASeqfinder
version 4.0.5
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

version="4.0.5"

### loading seqfinder and abricate results tables
if len(sys.argv)>1:
    df_abr = pd.read_csv(sys.argv[1], sep='\t')
    df_sf = pd.read_csv(sys.argv[2]) 
    df_original = pd.read_csv(sys.argv[3])
    efsa_dict = pd.read_csv(sys.argv[4])
    first_line = sys.argv[5]
    abricate_version = sys.argv[6]
else:  # just for developing code
    pd.set_option('display.max_columns', None)
    df_abr = pd.read_csv('/home/nickduggett/seqfinder_testing/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20230526/LREC5845/LREC5845.abricate', sep='\t')
    df_sf = pd.read_csv('/home/nickduggett/seqfinder_testing/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20230526/LREC5845/LREC5845_CompareTo_AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6_good_snps.csv')
    df_original= pd.read_csv('/home/nickduggett/seqfinder_testing/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6/20230526/LREC5845/LREC5845_CompareTo_AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6.csv')
    efsa_dict=pd.read_csv('/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.5/EFSA_panel/EFSA_antimcriobial_panel_dictionary_20230526.csv')
    abricate_version="1.0.1"
    with open('/home/nickduggett/seqfinder_testing/LREC5845.fasta') as f:
        first_line = f.readline().strip()

# Determining the versions of the assemblers used
match = re.search(r'Unicycler.*', first_line)
if match:
    assembler_versions = match.group()
else:
    assembler_versions = "Unknown"

# Extract the first version number
match = re.search(r'v([\d.]+)', assembler_versions)
unicycler_version = match.group(1) if match else "Unknown"

# Extract the second version number
match = re.search(r'SPAdes_v([\d.]+)', assembler_versions)
spades_version = match.group(1) if match else "Unknown"

###grouping plasmids by contig and choosing a representative 
#Make list of genes from the SeqFinder output
df_abr["LOCATION"]=df_abr['START'].astype(str)+".."+df_abr['END'].astype(str)

plasmids_only=df_abr[(df_abr.GENE.str.contains("PlasR"))]
grouped_plasmids = plasmids_only.loc[plasmids_only.groupby("SEQUENCE")["%COVERAGE"].idxmax()]
#### removed indexes with plasmids and merging the rest with the representatives
drop_plasmids=df_abr[(~df_abr.GENE.str.contains("PlasR"))]
merge_grouped_plasmids=drop_plasmids,grouped_plasmids
df_abr=pd.concat(merge_grouped_plasmids)


#Circumvents the filtering rule that would exclude ant3
ant3_rule = df_abr[(df_abr.GENE.str.contains("ant3") & (df_abr['%COVERAGE'] >=80))]
#Drop hits for chromosomal genes from the abricate output, ones less than 90% coverage and identity, and those that are plasmid genes
filtering_rule = df_abr[(~df_abr.GENE.str.contains("Chr") &(df_abr['%COVERAGE'] >= 90) &(df_abr['%IDENTITY'] >= 90) )]
#Merge ant3 back into the dataframe after circumventing the filter
merge_ant3 = filtering_rule,ant3_rule
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

abricate_dictionary_identity = filter_for_sf_genes[['GENE','%IDENTITY']]
abricate_dictionary_identity = abricate_dictionary_identity.set_index('GENE').T.to_dict('list')

abricate_dictionary_location = filter_for_sf_genes[['GENE','LOCATION']]
abricate_dictionary_location = abricate_dictionary_location.set_index('GENE').T.to_dict('list')



'''
#Map the dictionaries onto the SeqFinder dataframe into columns listed below
df_sf['result_abr']= df_sf['gene'].map(abricate_dictionary_gene)
df_sf['contig_abr']= df_sf['id'].map(abricate_dictionary_sequence)
df_sf['coverage_abr']= df_sf['id'].map(abricate_dictionary_coverage)
df_sf['identity_abr']= df_sf['id'].map(abricate_dictionary_identity)
df_sf['location_abr']= df_sf['id'].map(abricate_dictionary_location)

#Clean up the abricate columns to only include the text and no brackets
df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[2:-2])
df_sf['contig_abr'] = df_sf['contig_abr'].map(lambda x: str(x)[2:-2])
df_sf['coverage_abr'] = df_sf['coverage_abr'].map(lambda x: str(x)[1:-1])
df_sf['identity_abr'] = df_sf['identity_abr'].map(lambda x: str(x)[1:-1])
df_sf['location_abr'] = df_sf['location_abr'].map(lambda x: str(x)[2:-2])
#Replace 'a' (a carry over from NaN) with blanks
df_sf=df_sf.replace({'a':''})
'''
def clean_df_sf(df_sf, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_coverage, abricate_dictionary_identity, abricate_dictionary_location):
    # Map the dictionaries onto the SeqFinder dataframe into columns listed below
    df_sf['result_abr'] = df_sf['gene'].map(abricate_dictionary_gene)
    df_sf['contig_abr'] = df_sf['id'].map(abricate_dictionary_sequence)
    df_sf['coverage_abr'] = df_sf['id'].map(abricate_dictionary_coverage)
    df_sf['identity_abr'] = df_sf['id'].map(abricate_dictionary_identity)
    df_sf['location_abr'] = df_sf['id'].map(abricate_dictionary_location)

    # Clean up the abricate columns to only include the text and no brackets
    df_sf['result_abr'] = df_sf['result_abr'].map(lambda x: str(x)[2:-2])
    df_sf['contig_abr'] = df_sf['contig_abr'].map(lambda x: str(x)[2:-2])
    df_sf['coverage_abr'] = df_sf['coverage_abr'].map(lambda x: str(x)[1:-1])
    df_sf['identity_abr'] = df_sf['identity_abr'].map(lambda x: str(x)[1:-1])
    df_sf['location_abr'] = df_sf['location_abr'].map(lambda x: str(x)[2:-2])

    # Replace 'a' (a carry over from NaN) with blanks
    df_sf = df_sf.replace({'a': ''})
    
    return df_sf
df_sf = clean_df_sf(df_sf, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_coverage, abricate_dictionary_identity, abricate_dictionary_location)


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
'''
df_sf['plasmid_abr']= df_sf['contig_abr'].map(abricate_dictionary_plasmid)
df_sf['plasmid_abr'] = df_sf['plasmid_abr'].map(lambda x: str(x)[2:-2])
'''

#Filtering some pesky genes
values_to_check = ["TEM", "SHV", "CTX", "OXA", "CMY","tet","aadA"]
if df_sf[df_sf["id"].str.contains('|'.join(values_to_check)) & (df_sf["id"] == df_sf["result_abr"])].shape[0] > 0:
#if "betaL-g0285_CMY-2" in df_sf["id"].values:
    # Create a boolean mask to identify rows where "CMY" is an exact match
    mask = df_sf["id"].str.contains('|'.join(values_to_check))
    #df_sf = df_sf.loc[mask & (df_sf["id"] == df_sf["result_abr"]) | ~mask]
    if mask.any() and (df_sf.loc[mask, "id"] != df_sf.loc[mask, "result_abr"]).any():
        df_sf = df_sf.loc[mask & (df_sf["id"] == df_sf["result_abr"]) | ~mask]
else:
    df_sf=df_sf

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
            new_row['class'] = matching_row['class'].iloc[0]
            for col, mapping in column_mapping.items():
                new_row[col] = df_original[df_original['id'] == string][mapping].iloc[0]
            new_rows.append(new_row)
        return new_rows
    # Apply the function to rows with multiple strings and missing ids
    new_rows = []
    for idx, row in df_sf.iterrows():
        if row['multiple_strings'] and any(id_str not in unique_ids for id_str in row['result_abr']):
            new_rows.extend(create_missing_rows(row))
    
    # Append new rows to the DataFrame
    #df_sf = df_sf.append(new_rows, ignore_index=True) #Original code
    new_rows_df = pd.DataFrame(new_rows)
    df_sf = pd.concat([df_sf, new_rows_df], ignore_index=True)
    
    # Drop the intermediate column
    df_sf.drop('multiple_strings', axis=1, inplace=True)
    # Apply the functions from before that iterate through the dataframe again with the updated information
    df_sf = clean_df_sf(df_sf, abricate_dictionary_gene, abricate_dictionary_sequence, abricate_dictionary_coverage, abricate_dictionary_identity, abricate_dictionary_location)
    df_sf = process_plasmid_abr(df_sf, abricate_dictionary_plasmid)
    # Re-apply the efsa dictionary for the new row
    efsa_as_dict = efsa_dict.set_index('id').T.to_dict('list')
    df_sf['antimicrobial']= df_sf['id'].map(efsa_as_dict)
    # Formatting the dataframe so it looks a bit nicer
    df_sf=df_sf.fillna('')
    df_sf['antimicrobial'] = df_sf['antimicrobial'].map(lambda x: str(x)[2:-2])
    
#Reformat the output columns in the order we desire
abricate_columns=(df_sf.loc[:,['plasmid_abr','coverage_abr','identity_abr','location_abr']])
df_sf=df_sf.iloc[:,:-4]
merge=df_sf,abricate_columns

df_sf=pd.concat([df_sf, abricate_columns], axis=1, sort=False)

#Put SeqFinder version into final output
df_sf['SeqFinder_version']=version

df_sf['Abricate_version']=abricate_version

df_sf['Unicycler_version']=unicycler_version

df_sf['SPAdes_version']=spades_version

#Move columns around to make things look nicer
cols = df_sf.pop('reference'), df_sf.pop('EFSA_dict')
df_sf['reference'], df_sf['EFSA_dict'] = cols

#Getting rid of any "NaN"s carried over
df_sf.replace("NaN", np.nan, inplace=True)
df_sf.fillna("", inplace=True)

#Printing for developing only
if len(sys.argv)>1:
    output_filename = sys.argv[2][:-4]+'_abricate_seqfinder.csv'
    df_sf.to_csv(output_filename,index=False)
    print('Done! Saved output as {}.'.format(output_filename))
else:
    print(df_sf)


###To automatically solve the CTX, TEM and CMY multiple occurance issue this can be turned on
#id_list = ['CTX-','TEM-','CMY-']
#filter_for_TEM = df_sf.loc[df_sf['gene'].isin(id_list)]
#grouped_df = filter_for_TEM.groupby(by="gene")['snps'].min()




