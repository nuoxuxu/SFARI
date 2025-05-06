library(dplyr)
library(purrr)
library(ggplot2)

source("src/utils.R")

novel_splice_site_ClinVar <- read_vcf("export/variant/novel_splice_site_ClinVar.vcf", info=c("CLNVC", "CLNSIG", "MC", "GENEINFO"))

novel_splice_site_ClinVar %>%
    mutate(MC = str_remove(MC, ",.*$")) %>% 
    mutate(MC = str_extract(MC, "(?<=\\|)[^|]+")) %>%
    mutate(CLNSIG = map_chr(CLNSIG, ~ .x[1])) %>% 
    ggplot(aes(x=MC)) +
    geom_bar() +
    labs(x="Molecular consequences", y="Count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1))

novel_splice_site_ClinVar %>% 
    mutate(
        n_elements = str_count(MC, ",") + 1
    )

annovar_orfanage_res <- read_tsv("export/variant/annovar_with_orfanage_gtf/.hg38_multianno.txt")

annovar_gencode_res <- read_tsv("export/variant/annovar_with_gencode_gtf/.hg38_multianno.txt")

annovar_orfanage_res %>% 
    ggplot(aes(x=Func.refGene)) +
    geom_bar()

annovar_gencode_res %>% 
    ggplot(aes(x=ExonicFunc.refGene)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle=90, hjust=1))