library(dplyr)
library(ggplot2)
library(rtracklayer)
library(patchwork)


annot_peptides_hybrid <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble()

known <- annot_peptides_hybrid %>%
    filter(novelty == "known") %>% # filter to known peptides
    group_by(detected) %>% # group by detected status
    summarise(len = n(), .groups = "drop") %>% # count rows per group
    mutate(
        percent = len / sum(len) * 100, # calculate percentage over all groups
        type = "known" # add literal column type
    ) %>% 
    arrange(desc(detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent)

novel <- annot_peptides_hybrid %>%
    filter(novelty == "novel") %>% # filter to novel peptides
    group_by(detected) %>% # group by detected status
    summarise(len = n(), .groups = "drop") %>% # count rows per group
    mutate(
        percent = len / sum(len) * 100, # calculate percentage over all groups
        type = "novel" # add literal column type
    ) %>% 
    arrange(desc(detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent)

new_labels <- c("False" = "Not detected", "True" = "Detected")

p1 <- known %>% 
    ggplot(aes(x="", y=percent, fill=detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    scale_fill_manual(values = c("True"="#00BFC4", "False"="#F8766D"), labels = new_labels) +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 4) +
    theme(
        legend.position = "none"
    )

p2 <- novel %>% 
    ggplot(aes(x="", y=percent, fill=detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    scale_fill_manual(values = c("True"="#00BFC4", "False"="#F8766D"), labels = new_labels) +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 4) +
    labs(
        fill = "Detected"
    ) +
    theme(
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold")
    )

p1 + p2
ggsave("figures/figure_2/peptide_pie.pdf", width = 4, height = 2.5)
