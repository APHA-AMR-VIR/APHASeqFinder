#!/usr/bin/env python
''' 
Filter SNPs called by the APHA SeqFinder pipeline based on the following criteria:
'% mapped-gaps/real' == 100
'goodSnps' <= 30
'#non' > 0 (if gene in drop_non_syn_0)

Usage: goodSNPs.py csv_file [-h <help>]
'''

# Authors : Nicholas Duggett <Nick.Duggett@apha.gsi.gov.uk>
#           Nana Mensah <nmensah93@gmail.com>
# Created : 14/12/17
import sys
import pandas as pd
import logging
import warnings
pd.set_option('display.max_columns', 100)
# Configure a logger to print script messages to stdout. 
# Allows info messages to easily be printed, silenced or re-routed to a file.
# To switch off all messages, edit the command below with the argument: level='ERROR'.
logging.basicConfig(level='INFO', stream=sys.stdout, format="goodSNPs.py: %(message)s")

def amr_dataframe():
    '''Read dataframe from input csv file passed as command line argument. If no argument given,
    print usage string and exit.'''
    try:
        # Note: Filename is stored in sys.argv[1]
        return pd.read_csv(sys.argv[1])
        #return pd.read_csv('~/3820_CompareTo_AMRDatabase_20180308_ND.csv')
    # IOError occurs when file is not found. IndexError occurs when no arguments given.
    # Print usage string and exit for either event.
    except (IOError, IndexError):
        logging.info(__doc__)
        sys.exit()
# Create the name of the file that is being processed as a variable so it can be
# split into the strain name and added to the dataframe        
file_name = sys.argv[1]
#file_name='3820_CompareTo_AMRDatabase_20180308_ND.csv'
strain_name = file_name.split('_')
  
def gene_code(id_string):
    '''Given an id string, return the gene prefix, defined as the first four characters after the 
    first underscore.'''
    underscore_index = id_string.find('_')
    return id_string[underscore_index+1:underscore_index+5]

def filter_by_gene(dataframe):
    '''Return the input dataframe with a newly added gene column and row entries filtered to the 
    minimum goodSNPs in each gene category'''
    # Add a column to dataframe with the genes in which each SNP was called. Gene names are defined 
    # as the first 4 characters after the first underscore in the 'id' column of each entry
    dataframe['gene'] = [ gene_code(gene) for gene in dataframe['id'] ]
    # Add a column to the dataframe with the name of the strain the genes are from
    dataframe['strain']=strain_name[0]
    idlist = dataframe['id']
    amr_class=list()
    for i in idlist:
        amr_class.append(i[:5])
    dataframe['class']=amr_class
    # Create a dataframe of the lowest 'goodSnps' in each gene category
    grouped_df = dataframe.groupby(by="gene")['% mapped-gaps/real'].max()

    # Limit the dataframe to the columns 'gene' and 'good'Snps'. Using an unnamed (lambda) function,
    # create a list of 'True' or 'False' for each entry if it hold the minimum goodSNPs for its
    # gene category
    filter_criteria = dataframe[
                                ['gene','% mapped-gaps/real']
                                ].apply(lambda x: (x[0],x[1]) in grouped_df.iteritems(), axis=1)

    # Filter the dataframe by the True/False list and return
    return dataframe[filter_criteria]

def drop_non_syn(dataframe, non_syn_filter):
    '''Remove entries from dataframe where the 'gene' is in the list 'non_syn_filter' and there are 0
    non-synonymous SNPS present'''
    # Define a regular expression for the genes in non_syn_filter. Here each gene in the list is
    # separated by the or operator '|', with any character allowed either side. Eg '.*(gyrA|parC).*'
    gene_regex = ".*(" + "|".join(non_syn_filter) + ").*"

    # Ignore the warning printed when using pd.Series.str.contains() below
    warnings.filterwarnings("ignore", 'This pattern has match groups')

    # Search 'gene' column of dataframe for entries matching gene_regex and 0 non-synonymous SNPs
    matching_indices = dataframe.index[
                                       (dataframe['gene'].str.contains(gene_regex)) &
                                       (dataframe['#non'] == 0)
                                       ]
    # Return dataframe with entries that matched criteria (matching_indices) removed
    if not matching_indices.empty:
        # Note: The '~' in the dataframe slice is the logical 'NOT' in pandas
        return dataframe[~dataframe.index.isin(matching_indices)]
    else:
        # Return original dataframe if no entries are to be dropped
        return dataframe

def main():
    ''''''
    # Read input csv file to dataframe
    data_raw = amr_dataframe()
    # This line makes floR circumvent the filtering rules below if floR maps to more than 99%
    floR_rule = data_raw[(data_raw['% mapped-gaps/real'] >= 99) & (data_raw['id'] == 'chlor-g1585_floR')]
    gyrA_rule = data_raw[(data_raw['% mapped-gaps/real'] >= 98) & (data_raw['goodSnps'] <= 100) & (data_raw['id'] == 'quino-g1763_gyrA_ecol_Chr')]
    logging.info('Analysing results from "{}" '.format(sys.argv[1]))

    # Filter dataframe by % real SNPs and goodSNP count
    data_low_snps = data_raw[
                        (data_raw['% mapped-gaps/real'] >= 95) & 
                        (data_raw['goodSnps'] <= 100)
                       ].copy()
    # We now drop floR from the dataframe to ensure that it does not occur twice if the %mapped threshold is set below 99
    drop_floR = data_low_snps[(data_low_snps.id != 'chlor-g1585_floR') & (data_low_snps.id != 'quino-g1763_gyrA_ecol_Chr')]
    #drop_gyrA = data_low_snps[data_low_snps.id != 'quino-g1763_gyrA_ecol_Chr']
    # Merge the floR_rule above with the filtering rules set out in data_low_snps
    merge_floR = drop_floR,floR_rule,gyrA_rule
    data_low_snps = pd.concat(merge_floR)
    list_of_genes=['betaL-g0197_ampC_Chr','macro-g1725_mdfA_Chr','quino-g1763_gyrA_ecol_Chr','quino-g1809_parC_ecol_Chr','nitro-g1757_nfsA_Chr','fosfo-g1604_glpT_Chr','fosfo-g1606_ptsI_Chr','fosfo-g1607_uhpA_Chr','colis-g2212_phoP_EC_MG1655_Chr','colis-g2215_acrR_EC_MG1655_Chr']
    normalisation_genes = data_low_snps[data_low_snps.id.isin(list_of_genes)]
    number_of_genes = normalisation_genes.shape[0]
    sum_of_normalisation_genes = normalisation_genes['meanCov'].sum()
    mean_of_normalisation_genes = (sum_of_normalisation_genes / number_of_genes)
    data_low_snps["normCov"]=data_low_snps['meanCov'] / mean_of_normalisation_genes
    # Select dataframe row entries with the lowest number of SNPs for each 'gene' category
    data_gene_filtered = filter_by_gene(data_low_snps).copy()
    #Create filter rule to group each set of genes by the lowest number of goodsnps
    data_gene_filtered_snps=data_gene_filtered.groupby("gene")["goodSnps"].min()
    filter_criteria_snps = data_gene_filtered[['gene','goodSnps']].apply(lambda x: (x[0],x[1]) in data_gene_filtered_snps.iteritems(), axis=1)
    #Apply the filtered results to the main dataframe
    data_gene_filtered=data_gene_filtered[filter_criteria_snps]
    logging.info('Found {}/{} entries matching filter criteria'.format(data_gene_filtered.shape[0], 
                                                                       data_raw.shape[0]))    

    # Create a list of AMR genes to pass to drop_non_syn(). These strings will be searched for
    # in the 'gene' column, and will only be reported if non-synonymous SNPS are present.
    # Note: Each string in the list is searched with a wildcard before and after e.g *gyrA*
    drop_non_syn_0 = ['gyr','par','pmr','P3','P4','P5','sox','folP','bas','pts','Pc_P','Pa_P','omp','nfs','pho','etk','acr','amp','uhp','glp','mdf','mur']

    # Remove entries with 0 non-synonymous SNPs for each gene in drop_non_syn_0
    data_output = drop_non_syn(data_gene_filtered, drop_non_syn_0).copy()
    logging.info('Searched for #non == 0 in {}. '.format(drop_non_syn_0) +
                 'Report will display {} entries.'.format(data_output.shape[0]))
###Filter tet34 unless it hits at more than 90% mapped-gaps/real
    tet34_filter = data_output[(data_output['id'] == 'tetra-g1901_tet34') & (data_output['% mapped-gaps/real'] >= 90)]
    #drop_tet34 = data_output[(data_output.id != 'tetra-g1901_tet34')]
###Filter ampP from data_output if it does not contain snps in position 110 or 120
    ampP_filter = (data_output[data_output['gene'].str.contains("ampP")])
    ampP_snps_stringent = ['110-C-T-','120-']
    ampP_filter_snps = ampP_filter['annot'].str.contains('|'.join(ampP_snps_stringent))
    df_ampP_filter_snps = ampP_filter[ampP_filter_snps]
###Setting up variable to remove ampP, gyrA and parC from original output
    #drop_ampP_gyrA_parC = data_output[(data_output.gene != 'ampP') & (data_output.gene != 'gyrA') & (data_output.gene != 'parC') & (data_output.id != 'tetra-g1901_tet34')& (data_output.gene != 'ant3')]
    drop_ampP_gyrA_parC = data_output[~data_output.id.str.contains("_Chr")&(data_output.id != 'tetra-g1901_tet34')& (data_output.gene != 'ant3')]   
    data_output_chromosomal=(data_output[data_output['id'].str.contains('_Chr')])
    filter_list = ['gyrA', 'parC']
    df_gyrA_and_parC = (data_output[data_output.gene.isin(filter_list)])
###Filter the gyrA and parC dataframe to show only those that contain snps in the positions listed below and assign this to a new dataframe
    gyrA_parC_snps = ['248-C-T', '259-G-A-','239-G-T-']
    gyrA_parC_cip = df_gyrA_and_parC['annot'].str.contains('|'.join(gyrA_parC_snps))
    df_gyrA_parC_cip = df_gyrA_and_parC[gyrA_parC_cip]
###Filter multiple ant3, keeping just the longest version if both are present
    ant3_filter = (data_output[data_output['gene'].str.contains("ant3")])
    ant3_filter_long = ant3_filter.loc[ant3_filter.groupby('gene')['mapped-len'].idxmax()]
###Delete the original ampP, gyrA, parC, tet34 and ant3 rows from the original dataframe and merge ampP, gyrA or parC that passed the filters in ampP_snps_stringent or gyrA_parC_snps
    merge = drop_ampP_gyrA_parC,tet34_filter,df_gyrA_parC_cip,df_ampP_filter_snps,ant3_filter_long
    filtered_gyrA_parC_ampP = pd.concat(merge)
###Using gene and antimicrobial EFSA dictions to make another column with the AMR conferred by the gene
    test_dict_df = pd.read_csv('~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.1/seqfinder/EFSA_antimcriobial_panel_dictionary_260918.csv')
    df_as_dict = test_dict_df.set_index('id').T.to_dict('list')
    filtered_gyrA_parC_ampP['antimicrobial']= filtered_gyrA_parC_ampP['id'].map(df_as_dict)
    # Set output filename and columns
    output_filename = sys.argv[1].split('_')[0] + '_good_snps.csv'
    output_filename_chromosomal = sys.argv[1].split('_')[0] + '_good_snps_only_chromosomal.csv'
    original_columns = list(data_raw) # Get columns of input csv
    original_columns.insert(1, 'gene') # Set 'gene' as second column of output
    original_columns.insert(2, 'class') # Set 'class' as the third column of output
    original_columns.insert(0, 'strain') # Set 'strain' as the fourth column
    original_columns.insert(3, 'antimicrobial') # 'antimicrobial as the fifth column'
    original_columns.insert(8, 'normCov') # 'normCov as the ninth column'
    # Assign a list of column headers to be output in the final csv. This list can be edited here to
    # refine the final csv e.g. output_columns = ['id', '% mapped-gaps/real', 'goodSnps', '#syn', '#non']
    output_columns = original_columns
    output_df = filtered_gyrA_parC_ampP[output_columns].round(2)
    #print(output_df)
    # Write dataframe to CSV with suffix '_good_snps.csv'
    output_df.to_csv(output_filename, index=False)
    data_output_chromosomal.to_csv(output_filename_chromosomal, index=False)
    logging.info('Done! Output written to {}.'.format(output_filename))
    logging.info('Chromosomal output written to {}.'.format(output_filename_chromosomal))
if __name__=='__main__':
    main()
