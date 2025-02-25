# This script reformats the output of IUPred2A to be compatible with IsoformAnalyzeR
from src.utils import read_gtf

orfanage_gtf = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")
pbids = orfanage_gtf.unique("transcript_id")["transcript_id"].to_list()

with open("export/iupred2a_result.txt", "r") as f:
    iupred2a_result = [line.strip() for line in f if not line.startswith("#")]

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

out = split_and_insert(iupred2a_result, pbids)
out = [line for line in out if r"*" not in line]

with open("export/iupred2a_processed_result.txt", "w") as f:
    f.write("\n".join(out))