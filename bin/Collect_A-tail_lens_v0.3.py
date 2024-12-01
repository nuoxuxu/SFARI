#!/usr/bin/env python

"""
This script will query a SMRTLink instance based on a provided set
of credentials, and pull a variety of run metrics for longitudinal
tracking. Notably, if an output file of the same name exists, it
will filter results by existing metrics to enable it to run quickly
as new runs are introduced over time. Users have the option to 
over-write an existing file, but the default is to rename exiting
files by date.

And example of the credentials file is below:

Inputs
USER=user
PASSWD=password
SMRTLINK=smrtlink.host.domain (or IP address)
PORT=8243

"""

import argparse
import pandas as pd

# Command line options
parser = argparse.ArgumentParser(description='Python script to do \
                                 things. Need to add description of \
                                 credentials file')
parser.add_argument("-r","--readstatsfile",
					type=str,
					help="Path to the 'isoseq collapse' read stats file")
parser.add_argument("-f","--flncreport",
					type=str,
					help="Path to the reports.csv file produced by \
                        isoseq refine")
parser.add_argument("-o","--output",
					type=str,
                    help="Path to output file.")
args=parser.parse_args()


read_stats_file = args.readstatsfile
flnc_report_file = args.flncreport
output = args.output

#####################################################################
##################### Function Definitions ##########################
#####################################################################

def parse_read_stats_file(read_stats_file):
    ## Reads the isoseq collapse read stats file
    ## groups the read names by the isoform they are associated with
    ## returns grouped table
    collapsed_reads = pd.read_csv(read_stats_file, delimiter='\t')
    collapsed_reads['id'] = collapsed_reads['id'].astype(str)
    return(collapsed_reads)

def parse_flnc_report(flnc_report_file):
    ## Reads the isoseq refine 
    flnc_report = pd.read_csv(flnc_report_file)
    flnc_report['id'] = flnc_report['id'].astype(str)
    Atail_dict = dict(zip(flnc_report['id'], flnc_report['polyAlen']))
    return(Atail_dict)

def filter_and_aggregate_reads(collapsed_reads, Atail_dict):
    filtered_reads = collapsed_reads[collapsed_reads["id"].isin(list(Atail_dict.keys()))]
    filtered_reads = filtered_reads.groupby('pbid')['id'].agg(list).reset_index()
    return(filtered_reads)


def convert_reads_to_Atail(read_list, Atail_dict):
    Atail = [Atail_dict[read] for read in read_list]
    return(Atail)

def make_Atail_ouput(collapsed_reads, Atail_dict):
    collapsed_reads.iloc[:, 1] = collapsed_reads.iloc[:, 1].apply(convert_reads_to_Atail, args=(Atail_dict,))
    collapsed_reads.columns.values[1] = "Atail_len"
    collapsed_reads["Atail_len"] = collapsed_reads["Atail_len"].apply(lambda x: ','.join(map(str, x)))
    return(collapsed_reads)


#####################################################################
##################### Processing ##########################
#####################################################################
## Step1: Read inputs
print('*' * 30)
print("reading input files")
print('*' * 30)
collapsed_reads = parse_read_stats_file(read_stats_file)
Atail_dict = parse_flnc_report(flnc_report_file)

## Step2: Filter reads by the dictionary
## Which represent per-sample A tails (or could be inclusive of multiple samples)
print('*' * 30)
print("Filtering reads to those in flncs report")
print('*' * 30)
filtered_reads = filter_and_aggregate_reads(collapsed_reads, Atail_dict)

## Step3:
print('*' * 30)
print("Associating poly-A lengths")
print('*' * 30)
isoform_Atail_df = make_Atail_ouput(filtered_reads, Atail_dict)

## Step 4. Make final output
print('*' * 30)
print("writing output to " + output)
print('*' * 30)
isoform_Atail_df.to_csv(output,sep='\t', index=False)