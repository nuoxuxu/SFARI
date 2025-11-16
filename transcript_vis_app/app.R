library(ggplot2)
library(ggtranscript)
library(dplyr)
library(tidyr)
library(readr)
library(shiny)
library(bslib)
library(patchwork)
library(rtracklayer)
library(arrow)
library(dplyr)

my_theme <- theme_bw() + theme(
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_blank(),
  strip.text.y = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.ticks.y = element_blank(),
  legend.position = "top",
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14)
)

theme_set(my_theme)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC"
)

time_point_labels <- c(
  "mean_iPSC" = "t00",
  "mean_NPC"  = "t04",
  "mean_CN"   = "t30"
)

classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

lr_log2_cpm <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
  dplyr::rename(
    NPC_1_2 = NPC_1_3,
    NPC_3_2 = NPC_3_3,
    CN_1_1 = CN_1_2,
    CN_1_2 = CN_1_3
  ) %>%
  mutate(
    across(
      where(is.numeric), 
      ~ log2((.x / sum(.x) * 1e6) + 1)
    )
  ) %>% 
  mutate(
    mean_iPSC = rowMeans(select(., starts_with("iPSC_")), na.rm = TRUE),
    mean_NPC  = rowMeans(select(., starts_with("NPC_")), na.rm = TRUE),
    mean_CN   = rowMeans(select(., starts_with("CN_")), na.rm = TRUE)
  ) %>% 
  select(isoform, mean_iPSC, mean_NPC, mean_CN) %>% 
  pivot_longer(
    cols = c("mean_iPSC", "mean_NPC", "mean_CN"),
    names_to = "time_point",
    values_to = "abundance"
  ) %>% 
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = "isoform"
  ) %>% 
  dplyr::filter(structural_category %in% c("full-splice_match", "novel_not_in_catalog", "incomplete-splice_match", "novel_in_catalog")) %>% 
  mutate(
    structural_category = structural_category_labels[structural_category], 
    time_point = time_point_labels[time_point]
  ) %>% 
  mutate(
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC"))
  )

orfanage_gtf <- rtracklayer::import("nextflow_results/V47/orfanage/orfanage.gtf") %>% as.data.frame()
full_gtf <- rtracklayer::import("nextflow_results/V47/final_transcripts.gtf") %>% as.data.frame()
tx_no_CDS <- setdiff(pull(distinct(full_gtf, transcript_id), transcript_id), pull(distinct(orfanage_gtf, transcript_id), transcript_id))
gtf <- bind_rows(orfanage_gtf, full_gtf %>% filter(transcript_id %in% tx_no_CDS))

gtf <- gtf %>%
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = c("transcript_id" = "isoform")
  ) %>% 
  filter(structural_category %in% c("full-splice_match", "novel_not_in_catalog", "incomplete-splice_match", "novel_in_catalog")) %>% 
  mutate(
    structural_category = structural_category_labels[structural_category]
  ) %>% 
  mutate(
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC"))
  )

plot_ggtranscript <- function(gene_of_interest) {
  exons <- gtf %>%
    filter(
      type == "exon",
      associated_gene == {{ gene_of_interest }}
    )
  CDS <- gtf %>%
    filter(
      type == "CDS",
      associated_gene == {{ gene_of_interest }}
    )
  exons %>%
    ggplot(
      aes(xstart = start, xend = end, y = transcript_id, strand = strand)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id"),
    ) +
    geom_range(
      fill = "white",
      height = 0.25
    ) +
    geom_range(
      data = CDS,
      aes(fill = structural_category)
    ) +
    facet_grid(structural_category ~ ., scales = "free", space = "free")
}

plot_abundance <- function(gene_of_interest) {
  lr_log2_cpm %>% 
    filter(associated_gene == gene_of_interest) %>%
    ggplot(aes(x = abundance, y = isoform, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("log2(CPM + 1)") +
    ylab(gene_of_interest) +
    theme(
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 14),
    ) +
    scale_fill_manual(breaks = c("t00", "t04", "t30"), values = c("#e0a19c", "#e14bd5", "#ff0011")) +
    facet_grid(structural_category ~ ., scales = "free", space = "free") 
}

plot_combined <- function(gene_of_interest) {
  ggtranscript_plot <- plot_ggtranscript(gene_of_interest)
  abundance <- plot_abundance(gene_of_interest)
  ggtranscript_plot | abundance + plot_layout(widths = c(2, 1))
}

get_n_transcripts <- function(gene_of_interest) {
  classification %>%
    filter(associated_gene == gene_of_interest) %>%
    nrow()
}

ui <- page_sidebar(
  sidebar = sidebar(
    title = "Comparison",
    selectizeInput(
      "select_genes",
      "Select genes from list below:",
      choices = choices <- unique(gtf$associated_gene)
    )
  ),
  plotOutput("combined_plot")
)

server <- function(input, output, session) {
  output$combined_plot <- renderPlot(
    {
      plot_combined(input$select_genes)
    },
    height = reactive({
      35 * sum(get_n_transcripts(input$select_genes)) + 100
    })
  )
}

shinyApp(ui, server)