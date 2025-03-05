#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
from src.utils import read_gtf
import subprocess
from pathlib import Path
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser(description='Filter the scans to detected peptides only, also generate corresponding Identification.csv file')
    parser.add_argument('--mzXML_file', action='store', dest='mzXML_file', type=str, required=True)
    parser.add_argument("--annot_peptides", action='store', dest='annot_peptides', type=str, required=True)
    params = parser.parse_args()

    annot_peptides = read_gtf(params.annot_peptides, attributes = ["gene_id", "detected", "novelty", "transcript_id"])
    mzXML_file = params.mzXML_file

    file_name = Path(mzXML_file).stem

    identifications = annot_peptides\
        .filter(
            pl.col("detected") == "True",
            pl.col("novelty") == "novel"
        )\
        .unique("gene_id")\
        .with_columns(
            file_name = pl.col("gene_id").str.split("_").map_elements(lambda s: s[0], return_dtype=pl.String)
        )\
        .filter(
            pl.col("file_name") == file_name
        )\
        .with_columns(
            pl.col("gene_id").str.split("_").map_elements(lambda s: s[1], return_dtype=pl.String).cast(pl.Int32).alias("scan_number"),
            pl.col("transcript_id").str.split(".").map_elements(lambda s: s[1], return_dtype=pl.String)
        )\
        .sort("scan_number")\
        .rename({"scan_number": "Scan", "transcript_id": "Sequence"})\
        .select(["Scan", "Sequence"])

    # Filter for the detected scans based on the scan number
        
    scan_num_list = identifications\
        .select("Scan")\
        .to_pandas()["Scan"].to_list()

    subprocess.run(
        f"""sed 's/"RAW"/"RAWData"/g' {mzXML_file} > {mzXML_file.replace(".mzXML", "_fixed.mzXML")}""",
        shell=True)

    subprocess.run(
        f"""apptainer exec -B /gpfs/fs0/project/s/shreejoy/nxu/SFARI,${{PWD}} "${{NXF_SINGULARITY_CACHEDIR}}/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif" \
wine msconvert {mzXML_file.replace(".mzXML", "_fixed.mzXML")} \
--filter "peakPicking vendor msLevel=2" \
--filter "scanNumber {" ".join(scan_num_list)}" \
--filter "msLevel 2" \
--mzML \
--outfile {file_name}_filtered.mzML
    """,
        shell=True,
        capture_output=True)

    # Get the precursor charge

    import xml.etree.ElementTree as ET
    tree = ET.parse("file1_filtered.mzML")
    root = tree.getroot()
    ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}

    # Look for cvParam elements with the specific accession:

    xpath_expr = './/mzml:cvParam[@accession="MS:1000041"]'
    Charge = []
    for param in root.findall(xpath_expr, ns):
        # Check that name="charge state" to avoid collisions with other param(s)
        if param.get('name') == 'charge state':
            value = param.get('value')
            Charge.append(value)

    def find_all(a_str, sub):
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1: return
            yield start
            start += len(sub) or 1 # use start += 1 to find overlapping matches

    def get_TMT_K_string(peptide_seq):
        idx_list = list(find_all(peptide_seq, "K"))
        if len(list(idx_list)) == 0:
            return ""
        else:
            return ";".join(["TMT11plex:" + str(i) for i in idx_list])+ ";"

    def get_Cysteines_string(peptide_seq):
        idx_list = list(find_all(peptide_seq, "C"))
        if len(list(idx_list)) == 0:
            return ""
        else:
            return ";".join(["add_C_cysteine:" + str(i) for i in idx_list]) + ";"

    def get_full_mod_string(peptide_seq):
        return get_TMT_K_string(peptide_seq) + get_Cysteines_string(peptide_seq) + "TMT11plex:0"
        
    pl.concat([identifications, pl.DataFrame({"Charge": Charge})], how = "horizontal")\
        .with_columns(
            Modification = pl.col("Sequence").map_elements(lambda x: get_full_mod_string(x), return_dtype=pl.String)
        )\
        .select(["Scan", "Sequence", "Charge", "Modification"])\
        .write_csv(f"{file_name}_identifications.csv", include_header=True)
    
if __name__ == "__main__":
    main()