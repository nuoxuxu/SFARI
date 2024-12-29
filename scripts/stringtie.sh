# Stringtie needs ts tag from minimap2 to be able to detect splice junctions. merged.bam does not contain ts tags
# minimap2 does not accept PacBio Bam format, convert BAM to fastq first
sbatch -J bamtofastq -t 0-2:0 -N 1 -n 1 -o slurm_logs/bamtofastq.out --wrap="/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env/bin/bedtools bamtofastq -i data/long_read/LUO26876.20240516/CN_1_2_150PM_CELL1/m84090_240409_212646_s4.fail_reads.bcM0004.bam -fq CN_1_2_150PM_CELL1.fq"

sbatch -J stringtie -t 0-5:0 -N 1 -n 1 -o slurm_logs/stringtie.out --wrap="/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env/bin/stringtie -L proc/merged.bam -o stie.gtf -p 40"