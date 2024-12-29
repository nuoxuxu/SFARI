#!/bin/bash
#SBATCH --job-name=pfam_scan
#SBATCH --output=slurm_logs/pfam_scan.out
#SBATCH --time=1-0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba activate pfam_scan
export PERL5LIB=/home/s/shreejoy/nxu/tools/PfamScan:$PERL5LIB
/home/s/shreejoy/nxu/tools/PfamScan/pfam_scan.pl -fasta full_AA.fasta -dir /home/s/shreejoy/nxu/tools/Pfam_flat_files -cpu 40 -outfile full_pfam_scan_results.txt