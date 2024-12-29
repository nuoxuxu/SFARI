#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
from src.single_cell import SingleCell
import argparse

def main():
    parser = argparse.ArgumentParser(description='Get the list pbids that are in the high confidence set')
    parser.add_argument('--h5ad_file', action='store', dest='h5ad_file', type=str, required=True)
    parser.add_argument('--output', action='store', dest='output', type=str, required=True)
    params = parser.parse_args()

    lr_bulk = SingleCell(params.h5ad_file)
    lr_bulk.var["isoform"].to_frame().write_csv(params.output, include_header=False)

if __name__ == "__main__":
    main()    