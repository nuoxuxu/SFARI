library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(RColorBrewer)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
    )

theme_set(my_theme)

cell_summary <- read_csv("data/cell_Diff_summary.csv")

cell_summary %>%
    pivot_longer(c(OCT4, CTIP2), names_to = "marker", values_to = "percentage") %>% 
    ggplot(aes(x=time_point, y=percentage, fill=marker)) +
    geom_col(position="dodge2") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "black", size = 4) +
    scale_fill_brewer(palette = "Set2")

ggsave("figures/figure_1/cell_quantification.pdf", width = 3.8, height = 2.5)  
