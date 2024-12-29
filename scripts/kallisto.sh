#!/bin/bash
mamba activate lrp

kallisto index --index="GENCODE_v39" "${GENOMIC_DATA_DIR}GENCODE/gencode.v39.transcripts.fa"