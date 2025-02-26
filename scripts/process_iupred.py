# This script reformats the output of IUPred2A to be compatible with IsoformAnalyzeR
from src.utils import read_gtf
from io import StringIO
import polars as pl

def split_and_insert(list_A, list_B):
    result = []
    header_index = 0
    if header_index < len(list_B):
        result.append(">" + list_B[header_index])
        header_index += 1
    for item in list_A:
        result.append(item)
        if "*" in item and header_index < len(list_B):
            result.append(">" + list_B[header_index])
            header_index += 1
    return result

def split_and_read(list_a, marker="*"):
    blocks = []
    current_block = []
    
    for line in list_a:
        # If the line contains the marker, end the current block (if not empty)
        if marker in line:
            if current_block:
                blocks.append(current_block)
                current_block = []
        else:
            current_block.append(line)
    # Append the last block if it exists
    if current_block:
        blocks.append(current_block)
    
    # Convert each block (a list of strings) into a Polars DataFrame
    dfs = []
    for block in blocks:
        # Join the lines with newline characters
        block_str = "\n".join(block)
        # Use StringIO so that read_csv can read from the string as if it were a file
        df = pl.read_csv(StringIO(block_str), separator="\t", has_header = False, new_columns = ["POS", "RES", "IUPRED2", "ANCHOR2"])\
            .drop("POS")\
            .with_row_index(name = "POS")\
            .with_columns(
                POS = pl.col("POS") + 1
            )
        b = pl.DataFrame(
            {"POS":0,
            "RES":"*",
            "IUPRED2":0.0,
            "ANCHOR2":0.0},
            )\
            .with_columns(
                pl.col("POS").cast(pl.UInt32)
            )
        df = pl.concat([df, b], how = "vertical")  
        dfs.append(df)
    
    return dfs

# Get list of IUPred2A results

with open("export/iupred2a_result.txt", "r") as f:
    iupred2a_result = [line.strip() for line in f if not line.startswith("#")]

iupred2a_dfs = split_and_read(iupred2a_result)

temp = pl.concat(iupred2a_dfs, how = "vertical")\
    .with_columns(
        POS = pl.col("POS").cast(pl.String),
        RES = pl.col("RES").cast(pl.String),
        IUPRED2 = pl.col("IUPRED2").cast(pl.String),
        ANCHOR2 = pl.col("ANCHOR2").cast(pl.String)
    )\
    .to_pandas()\
    .to_numpy()

list_A = ["\t".join(row) for row in temp]

# Get header lines

pbids = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
    .unique("transcript_id")["transcript_id"]\
    .to_list()

out = split_and_insert(list_A, pbids)
out = [line for line in out if r"*" not in line]

with open("export/iupred2a_processed_result.txt", "w") as f:
    f.write("\n".join(out))