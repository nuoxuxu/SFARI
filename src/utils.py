import os
import polars as pl
import sys
import glob

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