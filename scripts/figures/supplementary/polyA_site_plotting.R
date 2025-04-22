library(readr)
library(dplyr)
library(ggplot2)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

structural_category_polyA_site_percentage <- read_csv("data/structural_category_polyA_site_percentage.csv")

structural_category_polyA_site_percentage %>% 
    filter(structural_category2 != "incomplete-splice_match") %>% 
    ggplot(aes(x = structural_category2, y = polyA_site_percentage)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Structural Category", y = "% within polyA site") +
    coord_flip()

ggsave("figures/supplementary/polyA_site.pdf", width = 4, height = 3)