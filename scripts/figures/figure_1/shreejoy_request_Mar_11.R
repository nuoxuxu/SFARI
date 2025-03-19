library(ggplot2)
library(dplyr)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
    )

theme_set(my_theme)

colorVector <- c(
    "Unchanged" = "#25222286",
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00"
)

df <- tibble(
    structural_category = c("Unchanged", "FSM", "NIC", "NNC"),
    n_DTE = c(57631, 36864, 33048, 21826)
) %>% 
    mutate(structural_category = factor(structural_category, levels = c("Unchanged", "FSM", "NIC", "NNC")))

df %>% 
    ggplot(aes(structural_category, n_DTE, fill = structural_category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("Structural Category", values = colorVector) +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(
        x = "",
        y = expression("DTEs (x" ~ 10^3 * ")")
        )

ggsave("figures/figure_1/DTE_structural_category.pdf", width = 4, height = 3)