import os
import polars as pl
import sys
import glob
from collections import defaultdict
from bx.intervals.intersection import Interval,IntervalTree

pl.enable_string_cache()

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")

def read_gff(file, attributes=["ID"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf"{attribute}=([^;]+)").alias(attribute) for attribute in attributes]
                )
    return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqid", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])\
        .with_columns(
            [pl.col("attributes").str.extract(rf"{attribute}=([^;]+)").alias(attribute) for attribute in attributes]
            ).drop("attributes")

def read_outfmt(file):
    return pl.read_csv(file, separator="\t", new_columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gstart", "gend", "sstart", "send", "evalue", "bitscore"])

def read_fasta(fasta_file, gencode = False):
    """
    Reads a FASTA file and converts it into a Polars DataFrame.
    Removes '*' from sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        polars.DataFrame: A DataFrame with 'transcript_id' and 'seq' columns.
    """
    sequences = []
    transcript_id = None
    seq = []
    
    with open(fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if transcript_id is not None:  # Save previous entry
                    sequences.append((transcript_id, "".join(seq).replace("*", "")))  # Strip '*'
                # Extract transcript_id (substring before the first space)
                if gencode is False:
                    transcript_id = line[1:].split(" ", 1)[0]
                else:
                    transcript_id = line[1:].split("|")[1]
                seq = []  # Reset sequence
            else:
                seq.append(line)  # Collect sequence lines

    # Add the last sequence
    if transcript_id is not None:
        sequences.append((transcript_id, "".join(seq).replace("*", "")))  # Strip '*'
    
    # Convert to Polars DataFrame
    df = pl.DataFrame(sequences, schema=["transcript_id", "seq"], orient="row")
    return df

def read_SJ(file):
    return pl.read_csv(
            file, separator="\t", 
            new_columns=["chrom", "start", "end", "strand", "motif", "annotated", "unique_reads", "multi_reads", "max_overhang"], 
            schema_overrides={"chrom": pl.String}
        )

def gtf_to_SJ(gtf):
    import numpy as np
    return gtf\
        .filter(
            pl.col("feature")!="transcript"
        )\
        .group_by("transcript_id", maintain_order=True)\
        .agg(
            pl.col("strand").unique().map_elements(lambda l: l[0], return_dtype=pl.String),
            pl.col("seqname").unique().map_elements(lambda l: l[0], return_dtype=pl.String),
            pl.col("start").map_elements(lambda l: np.sort(np.array(l)-1)[1:].tolist(), return_dtype=pl.List(pl.Int64)).alias("end"),
            pl.col("end").map_elements(lambda l: np.sort(np.array(l)+1)[:-1].tolist(), return_dtype=pl.List(pl.Int64)).alias("start")
        )\
        .explode("start", "end")\
        .rename({"seqname": "chrom"})\
        .filter(pl.col("start").is_null().not_())
        
def write_fasta(df: pl.DataFrame, id_col: str, seq_col: str, output_file: str):
    """
    Writes a FASTA file from a Polars DataFrame.
    
    Parameters:
        df (pl.DataFrame): Polars DataFrame containing sequence data.
        id_col (str): Column name for sequence IDs.
        seq_col (str): Column name for amino acid sequences.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for row in df.iter_rows(named=True):
            f.write(f">{row[id_col]}\n{row[seq_col]}\n")

def collapse_isoforms_to_proteoforms(gtf):
    return gtf\
        .filter(pl.col("feature")=="CDS")\
        .group_by("transcript_id")\
        .agg(
            pl.col("start"),
            pl.col("end")
        )\
        .group_by(["start", "end"])\
        .agg(
            pl.col("transcript_id")
        )\
        .with_columns(
            base_isoform = pl.col("transcript_id").map_elements(lambda x: x[0], return_dtype=pl.String)
        )\
        .explode("transcript_id")\
        .rename({"transcript_id": "isoform"})\
        .select("isoform", "base_isoform")

def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    strandedness="auto",
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    single_end=False,
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
    recursive=False,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[:sanitise_name_index]
            )
        return sample

    def get_fastqs(extension, recursive=False):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        search_path = f"*{extension}"
        if recursive:
            search_path = f"**/*{extension}"
        return sorted(glob.glob(os.path.join(fastq_dir, search_path), recursive=recursive))

    read_dict = {}

    ## Get read 1 files
    for read1_file in get_fastqs(read1_extension, recursive):
        sample = sanitize_sample(read1_file, read1_extension)
        if sample not in read_dict:
            read_dict[sample] = {"R1": [], "R2": []}
        read_dict[sample]["R1"].append(read1_file)

    ## Get read 2 files
    if not single_end:
        for read2_file in get_fastqs(read2_extension, recursive):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fastq_1", "fastq_2", "strandedness"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, read_1 in enumerate(reads["R1"]):
                    read_2 = ""
                    if idx < len(reads["R2"]):
                        read_2 = reads["R2"][idx]
                    sample_info = ",".join([sample, read_1, read_2, strandedness])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)    

class CAGEPeak:
    """
    A class to represent and query CAGE (Cap Analysis of Gene Expression) peaks from a BED file.

    Attributes
    ----------
    cage_bed_filename : str
        The filename of the BED file containing CAGE peak data.
    cage_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects containing intervals of peaks.

    Methods
    -------
    __init__(cage_bed_filename):
        Initializes the CAGEPeak object with the given BED filename and reads the BED file to populate the peaks.
    read_bed():
        Reads the BED file and populates the cage_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=10000):
        Queries the CAGE peaks to determine if a given position falls within a peak and calculates the distance to the nearest TSS.
    """
    def __init__(self, cage_bed_filename):
        self._validate_input(cage_bed_filename)
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def _validate_input(self, cage_bed_filename):
        if not cage_bed_filename.endswith('.bed'):
            raise ValueError("CAGE peak file must be in BED format.")
        if not os.path.exists(cage_bed_filename):
            raise FileNotFoundError("CAGE peak file does not exist.")
        
    def read_bed(self):
        for line in open(self.cage_bed_filename):
            raw = line.strip().split('\t')
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            tss0 = int((start0+end1)/2)
            self.cage_peaks[(chrom,strand)].insert(start0, end1, (tss0, start0, end1))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within a cage peak>, <nearest distance to TSS>
        If the distance is negative, the query is upstream of the TSS.
        If the query is outside of the peak upstream of it, the distance is NA
        """
        within_peak, dist_peak = 'FALSE', float('inf')
        peaks = self.cage_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for tss0, start0, end1 in peaks:
            # Checks if the TSS is upstream of a peak
            if (strand == '+' and start0 > query and end1 > query) or \
            (strand == '-' and start0 < query and end1 < query):
                continue
            # Checks if the query is within the peak and the distance to the TSS
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            d = (tss0 - query) * (-1 if strand == '-' else 1)
            w = 'TRUE' if within_out else 'FALSE'
            
            if not within_peak=='TRUE':
                within_peak, d = w, (tss0 - query) * (-1 if strand=='-' else +1)
                if within_peak == 'TRUE' or abs(d) < abs(dist_peak):
                    dist_peak = d
                
            else:
                d = (tss0 - query) * (-1 if strand=='-' else +1)
                if abs(d) < abs(dist_peak) and not(w == 'FALSE' and within_peak == 'TRUE'):
                    within_peak, dist_peak = w, d 
        if dist_peak == float('inf'):
            dist_peak = 'NA'
        return within_peak, dist_peak

def read_CAGE_peaks(CAGE_peak):
    if CAGE_peak:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        return CAGEPeak(CAGE_peak)
    return None

class PolyAPeak:
    """
    A class to represent and query polyA peaks from a BED file.

    Attributes
    ----------
    polya_bed_filename : str
        The filename of the BED file containing polyA peak information.
    polya_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects representing intervals of peaks.

    Methods
    -------
    __init__(polya_bed_filename)
        Initializes the PolyAPeak object with the given BED filename and reads the BED file to populate polyA peaks.
    read_bed()
        Reads the BED file and populates the polya_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=100)
        Queries the polyA peaks to determine if a given position falls within a specified search window of any peak.
    """
    def __init__(self, polya_bed_filename):
        self.polya_bed_filename = polya_bed_filename
        self.polya_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def read_bed(self):
        for line in open(self.polya_bed_filename):
            raw = line.strip().split()
            chrom = raw[0]
            start0 = int(raw[1])
            end1 = int(raw[2])
            strand = raw[5]
            self.polya_peaks[(chrom,strand)].insert(start0, end1, (start0, end1))

    def find(self, chrom, strand, query, search_window=100):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within some distance to polyA>, <distance to closest>
        - if downstream, + if upstream (watch for strand!!!)
        """
        assert strand in ('+', '-')
        within_polyA, dist_polyA = 'FALSE', 'NA'
        hits = self.polya_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for start0, end1 in hits:
            # Checks if the query is within the tail and the distance to the 5'
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            distance = start0 - query if strand == '+' else query - end1

            if within_out:
                within_polyA = 'TRUE'
            if dist_polyA == 'NA' or abs(distance) < abs(dist_polyA):
                dist_polyA = distance

        return within_polyA, dist_polyA

def read_polyA_peaks(polyA_peak):
    if polyA_peak:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        return PolyAPeak(polyA_peak)
    return None

