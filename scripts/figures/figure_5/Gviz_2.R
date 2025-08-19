library(Gviz)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(readr)
library(dplyr)
library(arrow)

# Read input files

orfanage_gtf <- import("nextflow_results/V47/orfanage/orfanage.gtf")
final_expression <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
    mutate(mean_expression = mean(c_across(where(is.numeric)), na.rm = TRUE))
AGO1 <- read_csv("data/katherine/AGO1_ex_Nuo.csv")

# Define functions

rep_tx_by_exp <- function(tx_name, exons) {
    rep_n <- final_expression %>%
        filter(isoform == tx_name) %>%
        pull(mean_expression) %>%
        round()
    GRangesList(rep(list(exons[[tx_name]]), rep_n)) %>% unlist()
}

rlelist_to_gr <- function(cov_rl) {
    stopifnot(is(cov_rl, "RleList"))
    out <- GRangesList(lapply(names(cov_rl), function(ch) {
        rl <- cov_rl[[ch]]
        ends <- cumsum(runLength(rl))
        starts <- c(1L, head(ends, -1L) + 1L)
        vals <- runValue(rl)

        keep <- vals != 0 # drop zero-coverage stretches
        if (!any(keep)) {
            return(GRanges())
        } # no signal on this chr

        GRanges(
            seqnames = ch,
            ranges = IRanges(start = starts[keep], end = ends[keep]),
            score = as.numeric(vals[keep])
        )
    }))
    unlist(out, use.names = FALSE)
}

getDataTrack <- function(tx_cluster_id) {
    tx_cluster <- AGO1 %>% dplyr::filter(cluster_id == tx_cluster_id)

    exons <- subset(orfanage_gtf, (mcols(orfanage_gtf)$transcript_id %in% tx_cluster$pb_id) & (mcols(orfanage_gtf)$type == "exon"))
    exons <- split(exons, exons$transcript_id)
    grl <- names(exons) %>% sapply(rep_tx_by_exp, exons = exons)
    cov <- do.call(c, unname(grl)) %>% GenomicRanges::coverage()
    cov_gr <- rlelist_to_gr(cov)
    DataTrack(range = cov_gr, genome = "hg38", name = tx_cluster_id, type = "histogram")
}

# Start plotting

DataTracks <- sapply(AGO1 %>% distinct(cluster_id) %>% pull(cluster_id), getDataTrack)

AnnotTrack <- orfanage_gtf %>%
    subset(
        (mcols(orfanage_gtf)$transcript_id %in% pull(distinct(AGO1, pb_id), pb_id)) & (mcols(orfanage_gtf)$type == "exon")
    ) %>%
    AnnotationTrack(
        genome = "hg38",
        name = "Gene models",
        fill = "lightblue",
        col = "black",
        shape = "box",
        group = .$transcript_id,
        groupAnnotation = "group"
    )

plotTracks(c(GenomeAxisTrack(), DataTracks, AnnotTrack))