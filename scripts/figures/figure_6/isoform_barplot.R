library(dplyr)
library(tidyr)
library(ggplot2)

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, vjust = 0.5, angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20)
    )

theme_set(my_theme)

isoformswitch <- readRDS("export/IsoformSwitchAnalyzeR/isoformswitch.rds")
isoformFeatures <- isoformswitch$isoformFeatures

isoformFeatures %>% 
  filter(
    isoform_id == novel_pb_id,
    condition_1 == "t00",
    condition_2 == "t04"
  ) %>% 
  pivot_longer(
    cols = c("iso_value_1", "iso_value_2"),
    names_to = "condition",
    values_to = "iso_value"
  ) %>%
  mutate(
    condition = recode(condition, "iso_value_1" = "t00", "iso_value_2" = "t30")
  ) %>% 
  ggplot(aes(x=condition, y=iso_value)) +
  geom_col() +
  labs(
    y = "Isoform Expression"
  ) +
  theme(
    axis.title.x = element_text(size = 25, color = "black")
  )

ggsave("figures/figure_6/isoform_barplot.pdf", width = 10, height = 8)  
