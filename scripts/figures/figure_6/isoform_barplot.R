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

isoformswitch <- readRDS("export/IsoformSwitchAnalyzeR/isoformswitch_part2_v2.rds")
isoformFeatures <- isoformswitch$isoformFeatures

novel_pb_id <- "PB.104608.73"
isoformFeatures %>% filter(isoform_id == novel_pb_id)

isoformFeatures %>% 
  filter(
    isoform_id == novel_pb_id,
    condition_1 == "t00",
    condition_2 == "t30"
  ) %>% 
  pivot_longer(
    cols = starts_with("iso_value_") | starts_with("iso_stderr_"),
    names_to = c(".value", "condition"),
    names_pattern = "iso_(value|stderr)_(\\d+)"
  ) %>%
  mutate(
    condition = recode(condition, "1" = "t00", "2" = "t30")
  ) %>% 
  ggplot(aes(x=condition, y=value)) +
  geom_col() +
  geom_pointrange(
    aes(ymin = value - stderr, ymax = value + stderr)
  )  +
  labs(y = "Isoform Fraction (IF)") +
  theme(
    axis.title.x = element_text(size = 25, color = "black")
  )

ggsave("figures/figure_6/isoform_barplot.pdf", width = 10, height = 8)  
