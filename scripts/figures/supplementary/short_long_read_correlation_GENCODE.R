txdb <- makeTxDbFromGFF("/project/s/shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf")

tx2gene <- select(txdb, keys = transcripts(txdb)$tx_name, columns = "GENEID", keytype = "TXNAME")

files <- Sys.glob("nextflow_results/salmon_GENCODE47/*/quant.sf")

names(files) <- str_split(files, "/") %>%
    map_chr(3)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM", txOut = TRUE, txIdCol = "isoform")

sr_expression <- txi$counts %>% as_tibble(rownames = "isoform")