import polars as pl
from src.utils import read_gtf
from src.ryp import r, to_py, to_r

r(
    """
library(GenomicRanges)
library(parallel)    
    """
)

transdecoder = pl.read_csv("full_nt.fasta.transdecoder.gff3", separator = "\t", new_columns=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])\
    .drop_nulls()\
    .filter(pl.col("type") == "CDS")\
    .rename({"seqid": "transcript_id"})

cds = read_gtf("proc/merged_collapsed.filtered.gtf", attributes=["transcript_id", "gene_id"]).filter(pl.col("feature") == "transcript")\
    .join(transdecoder["transcript_id", "start", "end"], on = "transcript_id", how = "left")\
    .with_columns(
        cds_start = pl.col("start") + pl.col("start_right") - 1,
        cds_end = pl.col("end") + pl.col("end_right") - 1
    )\
    .drop_nulls(["cds_start", "cds_end"])\
    .select(["seqname", "cds_start", "cds_end", "strand", "transcript_id"])\
    .rename({"cds_start": "start", "cds_end": "end"})
to_r(cds, "cds")

merged_collapsed_exons = read_gtf("proc/merged_collapsed.filtered.gtf", attributes=["transcript_id", "gene_id"]).filter(pl.col("feature") == "exon")\
    .select(["seqname", "start", "end", "strand", "transcript_id"])
to_r(merged_collapsed_exons, "merged_collapsed")

r(
    """
cds_gr <- GRanges(
    seqnames = cds$seqname,
    ranges = IRanges(start = cds$start, end = cds$end),
    strand = cds$strand,
    transcript_id = cds$transcript_id
)

merged_collapsed_gr <- GRanges(
    seqnames = merged_collapsed$seqname,
    ranges = IRanges(start = merged_collapsed$start, end = merged_collapsed$end),
    strand = merged_collapsed$strand,
    transcript_id = merged_collapsed$transcript_id
)

transcript_id_list <- mcols(cds_gr)$transcript_id

get_cds <- function(my_transcript_id) {
    gr <- intersect(cds_gr[mcols(cds_gr)$transcript_id==my_transcript_id], merged_collapsed_gr[mcols(merged_collapsed_gr)$transcript_id==my_transcript_id])
    mcols(gr)$transcript_id <- my_transcript_id
    gr
}

intersect_exon <- do.call(c, mclapply(transcript_id_list, get_cds, mc.cores=40))
saveRDS(intersect_exon, "proc/intersect_exon.rds")
    """
)