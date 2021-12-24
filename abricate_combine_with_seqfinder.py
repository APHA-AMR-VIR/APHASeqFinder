#!/usr/bin/env python2


'''
APHASeqfinder
version 4.0.0
submitted to github on 23/12/2021
Javier Nunez, AMR Team, Bacteriology (originally from Nicholas Duggett)
Animal and Plant Health Agency
This script takes abricate output, filters it, takes output, compares it seqfinder output and puts into a new spreadsheet
'''


# Authors : Nicholas Duggett <Nick.Duggett@apha.gov.uk>
# Created : 01/04/19
import sys
import pandas as pd
import warnings
#warnings.simplefilter(action='ignore', category=UserWarning)

### loading seqfinder and abricate results tables
df_abr = pd.read_csv(sys.argv[1], sep='\t')
df_sf = pd.read_csv(sys.argv[2])

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
#filtering_rule = df_abr[(~df_abr.GENE.str.contains("Chr") &(df_abr['%COVERAGE'] >= 90) &(df_abr['%IDENTITY'] >= 90) &(df_abr['DATABASE'] == 'AMRPlasmid'))]
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

#Map the dictionarys onto the SeqFinder dataframe into columns listed below
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

#Pick the plasmid rep gene with the highest coverage as the one to report
list_of_contigs = df_sf['contig_abr']
filter_for_abr_contigs = df_abr.loc[df_abr['SEQUENCE'].isin(list_of_contigs)]
filter_for_abr_contigs = filter_for_abr_contigs[(df_abr.GENE.str.contains("PlasR"))]
grouped_plasmids = df_abr.groupby(by="SEQUENCE")['%COVERAGE'].max()
abricate_dictionary_plasmid = filter_for_abr_contigs[['SEQUENCE','GENE']]
abricate_dictionary_plasmid = abricate_dictionary_plasmid.set_index('SEQUENCE').T.to_dict('list')


df_sf['plasmid_abr']= df_sf['contig_abr'].map(abricate_dictionary_plasmid)
df_sf['plasmid_abr'] = df_sf['plasmid_abr'].map(lambda x: str(x)[2:-2])

#Reformat the output columns in the order we desire
abricate_columns=(df_sf.loc[:,['plasmid_abr','coverage_abr','identity_abr','location_abr']])
df_sf=df_sf.iloc[:,:-4]
merge=df_sf,abricate_columns

df_sf=pd.concat([df_sf, abricate_columns], axis=1, sort=False)
output_filename = sys.argv[2][:-4]+'_abricate_seqfinder.csv'
df_sf.to_csv(output_filename,index=False)
print('Done! Saved output as {}.'.format(output_filename))


###To automatically solve the CTX, TEM and CMY multiple occurance issue this can be turned on
#id_list = ['CTX-','TEM-','CMY-']
#filter_for_TEM = df_sf.loc[df_sf['gene'].isin(id_list)]
#grouped_df = filter_for_TEM.groupby(by="gene")['snps'].min()




