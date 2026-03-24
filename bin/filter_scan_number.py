#!/usr/bin/env python3
import subprocess
from pathlib import Path
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser(description='Filter the scans to detected peptides only, also generate corresponding Identification.csv file')
    parser.add_argument('--mzXML_file', action='store', dest='mzXML_file', type=str, required=True)
    parser.add_argument("--novel_peptides", action='store', dest='novel_peptides', type=str, required=True)
    parser.add_argument("--percolator_res", action='store', dest='percolator_res', type=str, required=True)
    
    params = parser.parse_args()

    novel_peptides = pl.read_csv(params.novel_peptides)
    pep_to_scan_number = pl.read_csv(params.percolator_res, separator="\t")\
        .with_columns(
            pl.col("peptide").str.replace_all(r"M\[15.9949\]", "M")
        )\
        .with_columns(
            pep = pl.col("peptide").str.split(".").map_elements(lambda x: x[1], return_dtype=pl.String),
            prev_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[0], return_dtype=pl.String),
            next_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[2], return_dtype=pl.String)
        )\
        .filter(
            pl.col("q-value") < 0.05
        )\
        .unique("pep")\
        .select("PSMId", "pep")
    
    novel_peptides = novel_peptides\
        .join(
            pep_to_scan_number,
            left_on="original_pep",
            right_on="pep",
            how="left"
        )
    
    mzXML_file = params.mzXML_file

    file_name = Path(mzXML_file).stem

    identifications = novel_peptides\
        .unique("PSMId")\
        .with_columns(
            file_name = pl.col("PSMId").str.split("_").list.head(4).list.join("_")
        )\
        .filter(
            pl.col("file_name") == file_name
        )\
        .with_columns(Scan = pl.col("PSMId").str.split("_").list.get(4))\
        .rename({"original_pep": "Sequence"})\
        .select(["Scan", "Sequence"])

    # Filter for the detected scans based on the scan number
        
    scan_num_list = identifications\
        .select("Scan")\
        .to_pandas()["Scan"].to_list()

    subprocess.run(
        f"""sed 's/"RAW"/"RAWData"/g' {mzXML_file} > {mzXML_file.replace(".mzXML", "_fixed.mzXML")}""",
        shell=True)

    subprocess.run(
        f"""apptainer exec -B ${{PWD}} "${{NXF_SINGULARITY_CACHEDIR}}/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif" \
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
    tree = ET.parse(f"{file_name}_filtered.mzML")
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
            return ";".join(["add_K:" + str(i+1) for i in idx_list])+ ";"

    def get_Cysteines_string(peptide_seq):
        idx_list = list(find_all(peptide_seq, "C"))
        if len(list(idx_list)) == 0:
            return ""
        else:
            return ";".join(["add_C_cysteine:" + str(i+1) for i in idx_list]) + ";"

    def get_full_mod_string(peptide_seq):
        return get_TMT_K_string(peptide_seq) + get_Cysteines_string(peptide_seq) + "add_K:1"
        
    pl.concat([identifications, pl.DataFrame({"Charge": Charge})], how = "horizontal")\
        .with_columns(
            Modification = pl.col("Sequence").map_elements(lambda x: get_full_mod_string(x), return_dtype=pl.String)
        )\
        .select(["Scan", "Sequence", "Charge", "Modification"])\
        .write_csv(f"{file_name}_identifications.csv", include_header=True)
    
if __name__ == "__main__":
    main()

# ls *mzXML | xargs -I {} /scratch/nxu/SFARI/bin/filter_scan_number.py --mzXML_file {} --annot_peptides /scratch/nxu/SFARI/nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf