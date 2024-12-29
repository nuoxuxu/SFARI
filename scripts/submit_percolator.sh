#!/bin/bash
#SBATCH --job-name=submit_percolator
#SBATCH --output=slurm_logs/submit_percolator.out
#SBATCH --time=0-1:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

ls data/tc-1154/*.pin | xargs -I {} tail -n +2 {} > pooled.pin

echo "$(head -1 data/tc-1154/tc_1154_F1A_JM10561.pin)" | cat - pooled.pin > temp && mv -f temp pooled.pin

/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/percolator pooled.pin > pooled.tsv

awk '{
    for (i = 1; i <= NF; i++) {
        if (i <= 5) {
            printf "%s\t", $i;
        } else {
            printf "%s%s", $i, (i < NF ? "," : "");
        }
    }
    printf "\n";
}' OFS="\t" pooled.tsv > results/long_read/percolator.tsv