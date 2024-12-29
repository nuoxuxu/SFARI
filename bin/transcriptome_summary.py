#!/usr/bin/env python3
"""
This module prepares an isoform classification table and a gene info table 

Inputs:
--------------------------------------------------------------------------
 1. sqanti classification file 
 2. kalliso tpm file 
 3. kalliso ribodepletion tpm file 
 4. ensg_to_gene mapping file 
 5. enst_to_isoname mapping file 
 6. gene_len_stats table 
--------------------------------------------------------------------------

Outputs:
--------------------------------------------------------------------------
1. sqanti isoform table 
2. gene level info table 
--------------------------------------------------------------------------

"""

# Import Modules
import pandas as pd
import argparse

# Define Functions
def sqtab(sqanti_out, ensg_to_gene, enst_to_isoname):
    """
    Sorts data from Sqanti output 

    Note: We are not considering genes from sqanti output like ENSG00000242861.1_ENSG00000196187.12
    """

    # Import Data
    cols = ['isoform', 'length', 'structural_category','associated_gene','associated_transcript','subcategory']
    data = pd.read_csv(sqanti_out, delimiter="\t", usecols = cols)
    data.columns = ['pb_acc', 'len', 'cat', 'gene','transcript', 'cat2']

    # Map categories to acronyms and filter out anything that is not FSM, ISM, NNC or NIC
    data.replace({"novel_not_in_catalog":"NNC","novel_in_catalog":"NIC","incomplete-splice_match":"ISM","full-splice_match":"FSM"}, inplace = True)
    fdata = data[data.cat.isin(['FSM', 'ISM', 'NNC', "NIC"])]

    # Map gene -> gene_name
    gen_name = pd.read_csv(ensg_to_gene, delimiter="\t", header=None)
    gdict = pd.Series(gen_name.loc[:,1].values,index=gen_name.loc[:,0]).to_dict()
    df = fdata[['gene']]
    fdata['gene'] = fdata['gene'].map(gdict).fillna(df['gene'])

    # Drop cases like ENSG00000242861.1_ENSG00000196187.12 
    # fdata.drop(fdata[fdata['gene'] == df['gene']].index)

    # Map enst -> isoname
    trans = pd.read_csv(enst_to_isoname, delimiter="\t", header=None)
    tdict = pd.Series(trans.loc[:,1].values, index=trans.loc[:,0]).to_dict()
    df2 = fdata[['transcript']]
    fdata['transcript'] = fdata['transcript'].map(tdict).fillna(df2['transcript'])
    print("Isoform Table from sqanti output has been prepared")
    return fdata

def main():

    # Main Code 
    parser = argparse.ArgumentParser(description='Process transcriptome related input file locations')
    parser.add_argument('--sq_out', '-s', action='store', dest='sqanti_out', help = 'input : Sqanti Classification output location')
    parser.add_argument('--ensg_to_gene', '-gmap', action='store', dest='ensg_to_gene', help='ENSG -> Gene Map file location')
    parser.add_argument('--enst_to_isoname', '-imap', action='store', dest='enst_to_isoname', help='ENST -> Isoname Map file location')
    results = parser.parse_args()

    # Make Sqanti Isoform Table and output to a TSV
    sq_isotab = sqtab(results.sqanti_out, results.ensg_to_gene, results.enst_to_isoname)
    # sq_isotab['gene'] = sq_isotab['gene'].str.replace('_','-')

    # Make PB-Gene reference table
    pb_gene = sq_isotab[['pb_acc','gene']]
    # pb_gene.columns = ['isoform','gene']
    pb_gene = pb_gene.drop_duplicates()
    pb_gene.to_csv('pb_gene.tsv', sep="\t", index= False, na_rep='0')

if __name__ == "__main__":
    main()