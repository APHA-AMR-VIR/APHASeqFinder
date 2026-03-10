#!/usr/bin/env python

'''
APHASeqfinder
version 6.0.0
submitted to github on 23/12/2021
Javier Nunez, AMR Team, Bacteriology (originally from Nicholas Duggett)
Animal and Plant Health Agency
This script applies a series of filters to the seqfinder results
'''
#TODO: suggest moving this script to the end of the pipeline

## loading packages
import sys
import os
import logging
import warnings
import pandas as pd
import math

pd.set_option('display.max_columns', 100)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ANSI escape codes for text colors
class TextColours:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[39m'
    WARNING = '\033[93m'
    FAIL = '\033[38m'
    ENDC = '\033[0m'  # Reset color


## definition of functions
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
    dataframe['gene'] = [gene_code(gene) for gene in dataframe['id']]
    # Add a column to the dataframe with the name of the strain the genes are from
    strain=file_name.split(os.sep)[-2]
    dataframe['strain']=strain
    idlist = dataframe['id']
    amr_class=list()
    for i in idlist:
        amr_class.append(i[:5])
    dataframe['class']=amr_class
    # Create a dataframe of the lowest 'goodSnps' in each gene category
    grouped_df = dataframe.groupby(by="gene")['perc_mapped'].max()

    # Limit the dataframe to the columns 'gene' and 'good'Snps'. Using an unnamed (lambda) function,
    # create a list of 'True' or 'False' for each entry if it hold the minimum goodSNPs for its
    # gene category
    filter_criteria = dataframe[
                                ['gene','perc_mapped']
                                ].apply(lambda x: (x[0],x[1]) in grouped_df.items(), axis=1)

    # Filter the dataframe by the True/False list and return
    return dataframe[filter_criteria]

def drop_non_syn(dataframe, non_syn_filter):
    '''Remove entries from dataframe where the 'gene' is in the list 'non_syn_filter' and there are
    0 non-synonymous SNPS present'''
    # Define a regular expression for the genes in non_syn_filter. Here each gene in the list is
    # separated by the or operator '|', with any character allowed either side. Eg '.*(gyrA|parC).*'
    gene_regex = ".*(" + "|".join(non_syn_filter) + ").*"

    # Ignore the warning printed when using pd.Series.str.contains() below
    warnings.filterwarnings("ignore", 'This pattern has match groups')

    # Search 'gene' column of dataframe for entries matching gene_regex and 0 non-synonymous SNPs
    matching_indices = dataframe.index[
                                       (dataframe['gene'].str.contains(gene_regex)) &
                                       (dataframe['non'] == 0)
                                       ]
    # Return dataframe with entries that matched criteria (matching_indices) removed
    if not matching_indices.empty:
        # Note: The '~' in the dataframe slice is the logical 'NOT' in pandas
        return dataframe[~dataframe.index.isin(matching_indices)]
    else:
        # Return original dataframe if no entries are to be dropped
        return dataframe


### Error catching for genes not passing filter
class GenesFilterError(Exception):
    pass

######### Starting
######### reading arguments
if len(sys.argv)>1:
    file_name = sys.argv[1]
    per_ID=float(sys.argv[2])
    numsnps=float(sys.argv[3])
    efsa_dict=sys.argv[4]
    database_type=sys.argv[5]
    vir_dict=sys.argv[6]
    reference_name=sys.argv[7]
    bio_metal_dict=sys.argv[8]
    snp_pattern_dict=sys.argv[9]
else:  # just for developing code
    file_name = '/home/nickduggett/mnt/fsx-045/VM0533/VM0533J_SalAMR/furazolidone/test/AMR_with_plasmids_20250709/20250813/S02259-24/S02259-24_CompareTo_AMR_with_plasmids_20250709.csv'
    per_ID=70
    numsnps=100
    efsa_dict='/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_testing_area/APHASeqFinder_5.0.2/EFSA_panel/EFSA_antimcriobial_panel_dictionary_20230526.csv'
    database_type="AMR"
    vir_dict='/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_testing_area/APHASeqFinder_5.0.2/references/virulence/vir_dict_2022_06_17.csv'
    reference_name="AMR_with_plasmids_20230526"
    bio_metal_dict='/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_testing_area/APHASeqFinder_5.0.2/references/bio_metalinfectant/bio_metalinfectant_dictionary_2022_06_23.csv'
    snp_pattern_dict='/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_testing_area/APHASeqFinder_5.0.2/references/AMR/AMR_test_db_Chr_SNP_patterns.csv'
    
    
snp_ref = pd.read_csv(snp_pattern_dict).fillna("")
snp_ref['SNPs'] = snp_ref['snp'].str.split(",") + snp_ref['mut'].str.split(",")
snp_dict = dict(zip(snp_ref.gene, snp_ref.SNPs))
def snp_translator (gene,annot):
    if str(annot) == "nan":
        return("")
    annot = annot.replace(",", "")
    annot = annot.split(";")
    snplist = []
    for SNP in annot:
        if len(SNP) == 0:
            snplist.append("")
            break
        if "non" in SNP:
            wild = SNP.split(" ")[1].replace(']','').replace("'", "")
            mut = SNP.split(" ")[3].replace(']','').replace("'", "")
            pos = int(SNP.split("-")[0])
            pos = str(math.ceil(pos/3))
            #print("non")
            #print(wild+pos+mut)
            snplist.append(wild+pos+mut)
            '''
        else:
            wild = SNP.split("-")[1]
            mut = SNP.split("-")[2]
            pos = SNP.split("-")[0]
            #print("syn")
            #print("nt"+wild+pos+mut)
            #snplist.append("nt"+wild+pos+mut)
            snplist.append(wild+pos+mut)
            '''
    return(",".join(snplist))

def snp_patterner(gene, annot):
    foundlist = []
    if gene not in snp_dict:
        return "No SNPs listed for gene"
    if not isinstance(annot, str):
        return "No SNPs listed for gene"
    annot = annot.replace(",", "")
    annot = annot.split(";")
    snplist = []
    for SNP in annot:
        if len(SNP) == 0:
            snplist.append("")
            break
        if "non" in SNP:
            try:
                wild = SNP.split(" ")[1].replace(']', '').replace("'", "")
                mut = SNP.split(" ")[3].replace(']', '').replace("'", "")
                pos = int(SNP.split("-")[0])
                pos = str(math.ceil(pos / 3))
                snplist.append(wild + pos + mut)
            except (IndexError, ValueError):
                continue
        '''
        else:
            try:
                wild = SNP.split("-")[1]
                mut = SNP.split("-")[2]
                pos = SNP.split("-")[0]
                #snplist.append("nt" + wild + pos + mut)
                snplist.append(wild + pos + mut)
            except IndexError:
                continue
            '''
    for refsnp in snp_dict.get(gene, []):
        if refsnp in snplist:
            foundlist.append(refsnp)
    return ",".join(foundlist)


# Read input csv file to dataframe
try:
    data_raw=pd.read_csv(file_name)
    data_raw['reference'] = reference_name
except pd.errors.EmptyDataError:
    print("Warning: Empty CSV file "+file_name)
    file_name_failure=file_name+"failed"
    file_name_failure = file_name + "_failed.csv"
    with open(file_name_failure, "w", encoding="utf-8") as failure_file:
        failure_file.write("Empty CSV file: " + file_name)
    sys.exit() ####Might have to modify this
except (IOError, IndexError):
    print("Check "+file_name)
    with open(file_name_failure, "w", encoding="utf-8") as failure_file:
        failure_file.write("IOError or IndexError: " + file_name)
    sys.exit() ####Might have to modify this

if database_type=="AMR":
    data_raw["EFSA_dict"]=efsa_dict
elif database_type=="VIR":
    data_raw["Vir_dict"]=vir_dict
elif database_type=="BAM":
    data_raw["Bio_metal_dict"]=bio_metal_dict
else:
    pass #No action required if not using the three provided database types

if database_type=="AMR":
    #TODO: resolve these special cases generally
    # This line makes floR circumvent the filtering rules below if floR maps to more than 99% rather than the normal of 100 (this is legacy code now the threshold is set lower)
    floR_rule = data_raw[(data_raw['perc_mapped'] >= 98) & (data_raw['id'] == 'chlor-g1585_floR')]
    gyrA_rule = data_raw[(data_raw['perc_mapped'] >= 98) & (data_raw['good_snps'] <= 100) & (data_raw['id'] == 'quino-g1763_gyrA_ecol_Chr')]
    #logging.info('Analysing results from "{}" '.format(sys.argv[1]))

    # Filter dataframe by % real SNPs and goodSNP count
    data_low_snps = data_raw[(data_raw['perc_mapped'] >= per_ID) & (data_raw['good_snps'] <= numsnps)].copy()
    # We now drop floR from the dataframe to ensure that it does not occur twice if the %mapped threshold is set below 99
    drop_floR = data_low_snps[(data_low_snps.id != 'chlor-g1585_floR') & (data_low_snps.id != 'quino-g1763_gyrA_ecol_Chr')]
    #drop_gyrA = data_low_snps[data_low_snps.id != 'quino-g1763_gyrA_ecol_Chr']
    # Merge the floR_rule above with the filtering rules set out in data_low_snps
    merge_floR = drop_floR,floR_rule,gyrA_rule
    data_low_snps = pd.concat(merge_floR)
else:
    data_low_snps = data_raw[(data_raw['perc_mapped'] >= per_ID) & (data_raw['good_snps'] <= numsnps)].copy()

try:
    if len(data_low_snps) <1 or data_low_snps.empty:
        raise GenesFilterError(TextColours.WARNING+"No genes passed the filter set as "+str(per_ID)+"% mapped and "+str(numsnps)+" snps"+TextColours.ENDC)
except GenesFilterError as e:
    print("\n\nError: " + str(e) + "\n\nExiting\n\n")
    columns = ['warnings','strain','id','gene','antimicrobial','class','ref_len','mapped_len','mean_depth','norm_depth','non_calls','perc_mapped','snps','other','good_snps','syn','non','annotation','result_abr','contig_abr','coverage_abr','identity_abr','location_abr','plasmid_abr','coverage_map_abr','genes_on_same_contig','Platon_predicted_plasmid','Platon_RDS_likelihood','Generated_by','SeqFinder_version','Abricate_version','reference','EFSA_dict','MLST_reference']
    failed_isolate_df = pd.DataFrame(columns=dict.fromkeys(columns, []))
    failed_isolate_df.loc[0, 'strain'] = file_name.split(os.sep)[-2]
    failed_isolate_df = failed_isolate_df.fillna('No genes passed filter')
    output_filename = sys.argv[1].replace('.csv', '_good_snps.csv')
    failed_isolate_df.to_csv(output_filename, index=False)
    sys.exit()

# Select dataframe row entries with the lowest number of SNPs for each 'gene' category
data_gene_filtered = filter_by_gene(data_low_snps).copy()
'''
#Create filter rule to group each set of genes by the lowest number of goodsnps
data_gene_filtered_snps=data_gene_filtered.groupby("gene")["good_snps"].min()
#data_gene_filtered_snps=data_gene_filtered.groupby("gene")["non"].min()

filter_criteria_snps = data_gene_filtered[['gene','good_snps']].apply(lambda x: (x[0],x[1]) in data_gene_filtered_snps.items(), axis=1)
#Apply the filtered results to the main dataframe
data_gene_filtered=data_gene_filtered[filter_criteria_snps]

'''
# Step 1: Filter all rows with the minimum 'non' value per 'gene'
min_non_values = data_gene_filtered.groupby('gene')['non'].transform('min')
filtered_non = data_gene_filtered[data_gene_filtered['non'] == min_non_values]

# Step 2: From the above, filter all rows with the minimum 'syn' value per 'gene'
min_syn_values = filtered_non.groupby('gene')['syn'].transform('min')
filtered_syn = filtered_non[filtered_non['syn'] == min_syn_values]

# Final filtered DataFrame
data_gene_filtered = filtered_syn


data_output = data_gene_filtered.copy()
#TODO: WTF was I thinking, surely it'll be quicker to run snp_translator and then snp_patterner on that?
data_output['snps_causing_resistance'] = data_output.apply(lambda row: snp_patterner(row['id'],row['annotation']),axis=1)
data_output['annotation_non_synonymous_amino_acid'] = data_output.apply(lambda row: snp_translator(row['id'],row['annotation']),axis=1)
#data_output['better_name'] = data_output['id'].str.split("_").str[1] + "-" + data_output['snps_causing_resistance']
#data_output['better_name'] = data_output['better_name'].str.replace("No SNPs listed for gene","")
#data_output['better_name'] = data_output['better_name'].str.strip("-")

if database_type=="AMR":
    '''
    # Create a list of AMR genes to pass to drop_non_syn(). These strings will be searched for
    # in the 'gene' column, and will only be reported if non-synonymous SNPS are present.
    # Note: Each string in the list is searched with a wildcard before and after e.g *gyrA*
    drop_non_syn_0 = ['gyr','par','pmr','P3','P4','P5','sox','folP','bas','pts','Pc_P','Pa_P','omp','nfs','pho','etk','acr','amp','uhp','glp','mdf','mur']

    # Remove entries with 0 non-synonymous SNPs for each gene in drop_non_syn_0
    try:
        data_output = drop_non_syn(data_gene_filtered, drop_non_syn_0).copy()
    except AttributeError:
        strain=file_name.split(os.sep)[-2]
        print(TextColours.WARNING+"\n\nNo genes in your isolate ("+strain+") passed the filter set as "+str(per_ID)+"% mapped and "+str(numsnps)+" snps\n\n Exiting\n\n"+TextColours.ENDC)
        columns = ['warnings','strain', 'id', 'gene', 'antimicrobial',
                   'class', 'ref_len', 'mapped_len', 'mean_depth',
                   'norm_depth', 'non_calls', 'perc_mapped', 'snps',
                   'other', 'good_snps', 'syn', 'non', 'annotation',
                   'reference', 'EFSA_dict', 'result_abr', 'contig_abr',
                   'plasmid_abr', 'coverage_abr', 'identity_abr', 'location_abr']
        failed_isolate_df = pd.DataFrame(columns=dict.fromkeys(columns, []))
        failed_isolate_df.loc[0,'strain']=file_name.split(os.sep)[-2]
        failed_isolate_df = failed_isolate_df.fillna('No genes passed filter')
        output_filename = sys.argv[1].replace('.csv','_good_snps.csv')
        failed_isolate_df.to_csv(output_filename, index=False)
        sys.exit()
    logging.info('Searched for #non == 0 in {}. '.format(drop_non_syn_0) + 'Report will display {} entries.'.format(data_output.shape[0]))
    '''
    ###Filter tet34 unless it hits at more than 90% mapped-gaps/real
    #TODO: repalce special cases with general rules
    tet34_filter = data_output[(data_output['id'] == 'tetra-g1901_tet34') & (data_output['perc_mapped'] >= 90)]
    drop_tet34 = data_output[(data_output.id != 'tetra-g1901_tet34')]
    '''
    ###Filter ampP from data_output if it does not contain snps in position 110 or 120
    ampP_filter = data_output[data_output['gene'].str.contains("ampP")]
#This part is new
    ampP_filter['annotation'] = ampP_filter['annotation'].astype(str)
    ampP_snps_stringent = ['110-C-T-','120-']
    ampP_filter_snps = ampP_filter['annotation'].str.contains('|'.join(ampP_snps_stringent))
    df_ampP_filter_snps = ampP_filter[ampP_filter_snps]
    ###Setting up variable to remove ampP, gyrA and parC from original output
    #drop_ampP_gyrA_parC = data_output[(data_output.gene != 'ampP') & (data_output.gene != 'gyrA') & (data_output.gene != 'parC') & (data_output.id != 'tetra-g1901_tet34')& (data_output.gene != 'ant3')]
    drop_ampP_gyrA_parC = data_output[~data_output.id.str.contains("_Chr")&(data_output.id != 'tetra-g1901_tet34')& (data_output.gene != 'ant3')]
    '''
   
    data_output_chromosomal=(data_output[data_output['id'].str.contains('_Chr')])
    #output_filename ="good_snps_test.csv"
    #data_output_chromosomal.to_csv(output_filename, index=False)
    #print(data_output_chromosomal)
    data_output_chromosomal.loc[:, 'annotation'] = data_output_chromosomal['annotation'].fillna('')
    #Here we drop all ids that contain "_Chr" which should be chromosomal but we also need to get rid of ant3 because we have a special rule for that later
    drop_Chr = data_output[~data_output['id'].str.contains('_Chr') & (data_output.gene != 'ant3')]
    #print(drop_Chr)
    add_Chr = data_output_chromosomal[(data_output_chromosomal['snps_causing_resistance'] != "") & (data_output_chromosomal['snps_causing_resistance'] != "No SNPs listed for gene")] 
    #print(add_Chr)
    desired_order = ["strain", "id","gene", "class"] + [col for col in data_output_chromosomal.columns if col not in ["strain", "gene", "class"]]
    data_output_chromosomal=data_output_chromosomal[desired_order]
    
    '''
    # Create a new "warnings" column by applying the generate_warning function to each row
    filter_list = ['gyrA', 'parC']
    df_gyrA_and_parC = data_output[data_output.gene.isin(filter_list)]
    ###This part is new:
    df_gyrA_and_parC['annotation'] = df_gyrA_and_parC['annotation'].astype(str)
    ###Filter the gyrA and parC dataframe to show only those that contain snps in the positions listed below and assign this to a new dataframe
    gyrA_parC_snps = ['170-C-G','247-T-G', '248-C-A','248-C-T','259-G-A-','259-G-T','260-A-G-','239-G-T-','257-C-T','238-A-C-','250-G-A-']
    gyrA_parC_cip = df_gyrA_and_parC['annotation'].str.contains('|'.join(gyrA_parC_snps))
    df_gyrA_parC_cip = df_gyrA_and_parC[gyrA_parC_cip]
    '''
    ###Filter multiple ant3, keeping just the longest version if both are present
    ant3_filter = data_output[data_output['gene'].str.contains("ant3")]
    ant3_filter_long = ant3_filter.loc[ant3_filter.groupby('gene')['mapped_len'].idxmax()]
    ###Delete the original ampP, gyrA, parC, tet34 and ant3 rows from the original dataframe and merge ampP, gyrA or parC that passed the filters in ampP_snps_stringent or gyrA_parC_snps
    #merge = drop_ampP_gyrA_parC,tet34_filter,df_gyrA_parC_cip,df_ampP_filter_snps,ant3_filter_long
    merge = drop_Chr,add_Chr,tet34_filter,ant3_filter_long
    filtered_gyrA_parC_ampP = pd.concat(merge)
    #Rules to get rid of false positives that occur with some genes
    #TODO: resolve this problem generally with renamed gene groups for running purposes
    if 'betaL-g0444_CTX-M-55' in filtered_gyrA_parC_ampP['id'].values or 'betaL-g0362_CTX-M-1' in filtered_gyrA_parC_ampP['id'].values:
        mask = filtered_gyrA_parC_ampP['id'].str.contains('KLUC|KLUG')
        filtered_gyrA_parC_ampP.loc[mask & (filtered_gyrA_parC_ampP['non'] == '0'), 'non'] = '0'
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~(mask & (filtered_gyrA_parC_ampP['non'] != '0'))]
    else:
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~filtered_gyrA_parC_ampP['id'].str.contains('KLUC|KLUG')]
    if 'betaL-g0285_CMY-2' in filtered_gyrA_parC_ampP['id'].values:
        mask = filtered_gyrA_parC_ampP['id'].str.contains('BIL-1|LAT-1|CFE-1')
        filtered_gyrA_parC_ampP.loc[mask & (filtered_gyrA_parC_ampP['non'] == '0'), 'non'] = '0'
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~(mask & (filtered_gyrA_parC_ampP['non'] != '0'))]
    else:
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~filtered_gyrA_parC_ampP['id'].str.contains('BIL-1|LAT-1|CFE-1')]
    if 'betaL-g2446_SHV-12' in filtered_gyrA_parC_ampP['id'].values:
        mask = filtered_gyrA_parC_ampP['id'].str.contains('LEN-6|OHIO-1|OKP-A-11|OKP-A-12|OKP-A-13|OKP-A-15|OKP-B-15|OKP-B-17')
        filtered_gyrA_parC_ampP.loc[mask & (filtered_gyrA_parC_ampP['non'] == '0'), 'non'] = '0'
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~(mask & (filtered_gyrA_parC_ampP['non'] != '0'))]
    else:
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~filtered_gyrA_parC_ampP['id'].str.contains('LEN-6|OHIO-1|OKP-A-11|OKP-A-12|OKP-A-13|OKP-A-15|OKP-B-15|OKP-B-17')]
    if 'betaL-g0491_DHA-1' in filtered_gyrA_parC_ampP['id'].values:
        mask = filtered_gyrA_parC_ampP['id'].str.contains('MOR-2')
        filtered_gyrA_parC_ampP.loc[mask & (filtered_gyrA_parC_ampP['non'] == '0'), 'non'] = '0'
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~(mask & (filtered_gyrA_parC_ampP['non'] != '0'))]
    else:
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~filtered_gyrA_parC_ampP['id'].str.contains('MOR-2')]
        '''

        # Define rules to get rid of false positives that occur with some genes as a list of dictionaries
        rules = [
            {
                'trigger': ['betaL-g0444_CTX-M-55', 'betaL-g0362_CTX-M-1'],
                'filter': 'KLUC|KLUG'
            },
            {
                'trigger': ['betaL-g0285_CMY-2'],
                'filter': 'BIL-1|LAT-1|CFE-1'
            },
            {
                'trigger': ['betaL-g2446_SHV-12'],
                'filter': 'LEN-6|OHIO-1|OKP-A-11|OKP-A-12|OKP-A-13|OKP-A-15|OKP-B-15|OKP-B-17'
            },
            {
                'trigger': ['betaL-g0491_DHA-1'],
                'filter': 'MOR-2'
            }
        ]

        # Apply each rule
        for rule in rules:
            if any(gene in filtered_gyrA_parC_ampP['id'].values for gene in rule['trigger']):
                mask = filtered_gyrA_parC_ampP['id'].str.contains(rule['filter'])
                filtered_gyrA_parC_ampP.loc[mask & (filtered_gyrA_parC_ampP['non'] == '0'), 'non'] = '0'
                filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~(mask & (filtered_gyrA_parC_ampP['non'] != '0'))]
            else:
                filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~filtered_gyrA_parC_ampP['id'].str.contains(rule['filter'])]
        '''
    # Check if both "chlor-g1576_cml" and "chlor-g1578_cmlA1" are present in the "id" column
    if all(elem in filtered_gyrA_parC_ampP["id"].values for elem in ["chlor-g1576_cml", "chlor-g1578_cmlA1"]):
        # Filter the dataframe for rows containing either "chlor-g1576_cml" or "chlor-g1578_cmlA1"
        mask = filtered_gyrA_parC_ampP["id"].isin(["chlor-g1576_cml", "chlor-g1578_cmlA1"])
        filtered_df = filtered_gyrA_parC_ampP[mask]
        # Find the row with the highest "perc_mapped" value
        max_perc_mapped = filtered_df["perc_mapped"].max()
        # Check if both "chlor-g1576_cml" and "chlor-g1578_cmlA1" have the same highest "perc_mapped" value
        if (filtered_df["perc_mapped"] == max_perc_mapped).sum() > 1:
            # Remove "chlor-g1576_cml" and keep "chlor-g1578_cmlA1"
            filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP[~((filtered_gyrA_parC_ampP["id"] == "chlor-g1576_cml") & (filtered_gyrA_parC_ampP["perc_mapped"] == max_perc_mapped))]
        else:
            filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP
    else:
        filtered_gyrA_parC_ampP = filtered_gyrA_parC_ampP

else:
    filtered_gyrA_parC_ampP=data_gene_filtered

if database_type=="AMR":
    ###Using gene and antimicrobial EFSA dictions to make another column with the AMR conferred by the gene
    #First we take just the first two columns of the EFSA dict as that is all we want to include in the final df
    efsa_dict_df = pd.read_csv(efsa_dict, usecols=['id', 'antimicrobial'])
    #We convert the df into a dictionary
    df_as_dict = dict(zip(efsa_dict_df['id'], efsa_dict_df['antimicrobial']))
    #We finally apply the dictionary to the dataframe using the id value as a key
    filtered_gyrA_parC_ampP['antimicrobial']= filtered_gyrA_parC_ampP['id'].map(df_as_dict)
    # Set output columns
    output_columns = list(data_raw)
    output_columns.insert(1,'gene')
    output_columns.insert(2,'class')
    output_columns.insert(0,'strain')
    output_columns.insert(16,'snps_causing_resistance')
    output_columns.insert(17,'annotation_non_synonymous_amino_acid')
    output_columns.insert(3,'antimicrobial')
    #filtered_gyrA_parC_ampP['antimicrobial'] = filtered_gyrA_parC_ampP['antimicrobial'].map(lambda x: str(x)[2:-2])
    #Replace 'a' (a carry over from NaN) with blanks
    filtered_gyrA_parC_ampP=filtered_gyrA_parC_ampP.replace({'a':''})
elif database_type=="VIR":
    ###Using gene and Virulence dictionary to make another column with the Virulence function conferred by the gene
    vir_dict_df = pd.read_csv(vir_dict)
    columns=['Virulence_function','Associated_pathotypes']

    vir_dicts =  {col: vir_dict_df.set_index('id')[col].to_dict() for col in columns}
    for col in columns:
        filtered_gyrA_parC_ampP[col] = filtered_gyrA_parC_ampP['id'].map(vir_dicts[col]).map(str)
    '''
    ###Select just the columns "id" and "Virulence function" and set to new dataframe "vir_dict_gene_function"
    vir_dict_gene_function = vir_dict_df[['id','Virulence_function']]
    ###Change this dataframe to a dictionary list that can be called later
    vir_dict_gene_function = vir_dict_gene_function.set_index('id').T.to_dict('list')
    ###Select just the columns "id" and "Associated_pathotypes" and set to new dataframe "vir_dict_pathotypes"
    vir_dict_pathotypes = vir_dict_df[['id','Associated_pathotypes']]
    ###Change this datafram to a dictionary that can be called later
    vir_dict_pathotypes = vir_dict_pathotypes.set_index('id').T.to_dict('list')
    ###Make two new columns in our Seqfinder output dataframe called "Virulence_function" and "Associated_pathotypes"
    ###using the dictionaries we created to search for 'id' (the genes in our isolate)
    ###and match them with the associated keys in our dictionaries
    filtered_gyrA_parC_ampP['Virulence_function']= filtered_gyrA_parC_ampP['id'].map(vir_dict_gene_function)
    filtered_gyrA_parC_ampP['Associated_pathotypes']= filtered_gyrA_parC_ampP['id'].map(vir_dict_pathotypes)
    #Clean up the columns to only include the text and no brackets
    filtered_gyrA_parC_ampP['Virulence_function'] = filtered_gyrA_parC_ampP['Virulence_function'].map(lambda x: str(x)[2:-2])
    filtered_gyrA_parC_ampP['Associated_pathotypes'] = filtered_gyrA_parC_ampP['Associated_pathotypes'].map(lambda x: str(x)[2:-2])
    #Replace 'a' (a carry over from NaN) with blanks
    filtered_gyrA_parC_ampP=filtered_gyrA_parC_ampP.replace({'a':''})
    '''
    # Set output columns
    output_columns = list(data_raw)
    output_columns.insert(1,'gene')
    output_columns.insert(2,'Virulence_function')
    output_columns.insert(0,'strain')
    output_columns.insert(3,'Associated_pathotypes')
elif database_type == "BAM":
    bio_metal_dict_df = pd.read_csv(bio_metal_dict)
    columns = ['class', 'gene_function', 'sub_class', 'source_genus']
    
    bio_metal_dicts = {col: bio_metal_dict_df.set_index('id')[col].to_dict() for col in columns}

    for col in columns:
        filtered_gyrA_parC_ampP[col] = filtered_gyrA_parC_ampP['id'].map(bio_metal_dicts[col]).map(str)

    #Replace 'a' (a carry over from NaN) with blanks
    #filtered_gyrA_parC_ampP=filtered_gyrA_parC_ampP.replace({'a':''})
    # Set output columns
    output_columns = list(data_raw)
    output_columns.insert(1,'gene')
    output_columns.insert(2,'class')
    output_columns.insert(0,'strain')
    output_columns.insert(3,'sub_class')
    output_columns.insert(4,'gene_function')
    output_columns.insert(5,'source_genus')


else:
    output_columns = list(data_raw)
    output_columns.insert(1,'gene')
    output_columns.insert(0,'strain')
# Assign a list of column headers to be output in the final csv. This list can be edited here to
output_df = filtered_gyrA_parC_ampP[output_columns].round(2)

# filtering rows by the per_ID when database_type=="nonAMR"
if database_type!="AMR":
    output_df=output_df.loc[output_df['perc_mapped'] > per_ID]


###This part catches samples that have chromosomal genes above the threshold but do not have a known phenotype
try:
    if output_df.empty:
        raise GenesFilterError("No genes passed the filter set as "+str(per_ID)+"% mapped and "+str(numsnps)+" snps")
except GenesFilterError as e:
    print("\n\nError: " + str(e) + "\n\nExiting\n\n")
    columns = ['warnings','strain','id','gene','antimicrobial','class','ref_len','mapped_len','mean_depth','norm_depth','non_calls','perc_mapped','snps','other','good_snps','syn','non','annotation','result_abr','contig_abr','coverage_abr','identity_abr','location_abr','plasmid_abr','coverage_map_abr','genes_on_same_contig','Platon_predicted_plasmid','Platon_RDS_likelihood','Generated_by','SeqFinder_version','Abricate_version','reference','EFSA_dict','MLST_reference']
    failed_isolate_df = pd.DataFrame(columns=dict.fromkeys(columns, []))
    failed_isolate_df.loc[0, 'strain'] = file_name.split(os.sep)[-2]
    failed_isolate_df = failed_isolate_df.fillna('No genes passed filter')
    output_filename = sys.argv[1].replace('.csv', '_good_snps.csv')
    failed_isolate_df.to_csv(output_filename, index=False)
    sys.exit()

output_df["warnings"]=""
# Reorder the columns to make the "warnings" column the first column
output_df = output_df[["warnings"] + output_df.columns[:-1].tolist()]
output_df["annotation"] = output_df["annotation"].fillna("")

# Write dataframe to CSV with suffix '_good_snps.csv'
if len(sys.argv)>1:
    output_filename = sys.argv[1].replace('.csv','_good_snps.csv')
    output_df.to_csv(output_filename, index=False)
    if database_type=="AMR":
        output_filename_chromosomal =  sys.argv[1].replace('.csv','_good_snps_only_chromosomal.csv')
        data_output_chromosomal.to_csv(output_filename_chromosomal, index=False)
else:
    print(output_df)
    output_filename ="good_snps_test.csv"
    #output_df.to_csv(output_filename, index=False)
    print("Done!")
