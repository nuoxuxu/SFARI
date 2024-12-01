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
from pathlib import Path
import polars as pl

# Command line options
parser = argparse.ArgumentParser(description='Python script to do \
                                 things. Need to add description of \
                                 credentials file')
parser.add_argument("-r","--readstatsfile",
					type=str,
					help="Path to the 'isoseq collapse' read stats file")
parser.add_argument("-p","--prefix",
					type=str,
					help="Path to the directory that contains flnc.report.csv files")
parser.add_argument("-o","--output",
					type=str,
                    help="Path to output file.")
args=parser.parse_args()

#####################################################################
##################### Function Definitions ##########################
#####################################################################

def get_agg_polyALength(path):
    sample_name = str(path).split("/")[-3].rsplit("_", 2)[0]
    flnc_report = pl.read_csv(path)
    return flnc_report\
        .join(collapsed_reads, on='id', how='inner')\
        .group_by('pbid')\
        .agg(pl.col('polyAlen'))\
        .rename({"polyAlen": sample_name})

#####################################################################
##################### Processing ####################################
#####################################################################

collapsed_reads = pl.read_csv(args.readstatsfile, separator="\t")
path_list = list(Path(args.prefix).glob("*/outputs/flnc.report.csv"))

df = get_agg_polyALength(path_list[0])

for path in path_list[1:]:
    df = df.join(get_agg_polyALength(path), on='pbid', how='full', coalesce = True)

df.write_parquet(args.output)