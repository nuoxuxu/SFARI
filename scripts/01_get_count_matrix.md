Since QNAME of the input files to `isoseq collapse`, `mapped.bam` are designated by the sample name ("CN_2_2_transcript/4694701"), it is possible to merge all `mapped.bam` files for all samples and run `isoseq collapse` on this merged bam file. Memory should be ok because the size of a `mapped.bam` is less than 2 GB.

Once we get the mapping of the collapsed transcript names to the transcript names in each sample, we can get a count matrix of the collapsed transcripts in each of the samples.

Merge mapped.bam files, these bam files are outputs from `pbmm2`

```bash
sbatch  -t 0-5:0 --nodes=1 -J merge_bam -o merge_bam.out --wrap="ls data/long_read/LUO26876.20240514/*/outputs/mapped.bam | xargs samtools merge -o proc/merged.bam -@ 40"
```

Run isoseq collapse

```bash
sbatch -t 0-12:0 --nodes=1 -o isoseq_collapse.out -J isoseq_collapse --wrap="isoseq collapse -j 50 proc/merged.bam proc/merged_collapsed.gff"
```

Get mapping from ccs read name to sample name

```bash
for file in data/long_read/LUO26876.20240514/*/outputs/flnc.bam; do
	if [[ -f "$file" ]]; then
		first_field=$(samtools view "$file" | head -1 | cut -f1)
		if [[ -n "$first_field" ]]; then
			echo -e "$first_field\t$file" >> proc/id_to_sample.txt
		else
			echo -e "No data\t$file" >> proc/id_to_sample.txt
		fi
	else
		echo -e "File not found\t$file" >> proc/id_to_sample.txt
	fi
done
```

Filter out transcripts affected by intra-priming and RT switching

```bash
mamba activate tama
pigeon filter proc/merged_collapsed_classification.txt
```