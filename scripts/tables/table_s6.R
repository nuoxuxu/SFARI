library(dplyr)
library(readr)

All_variants <- read_tsv("export/variant/All_variants_used_in_project.tsv")
All_variants %>% 
    filter(variant_type == "denovo") %>% 
    select(-is_conserved) %>% 
    write_tsv("export/variant/Table_S6.tsv")