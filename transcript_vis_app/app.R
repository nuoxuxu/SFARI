library(ggplot2)
library(ggtranscript)
library(dplyr)
library(tidyr)
library(readr)
library(bslib)
library(patchwork)

combined_gtf <- read_csv("./shiny.csv")
pbid_abundance <- read_csv("./pbid_abundance.csv")
talon_abundance <- read_csv("./talon_abundance.csv")

ui <- page_fillable(
  selectizeInput(
    "select_genes",
    "Select genes from list below:",
    choices = choices <- unique(combined_gtf$gene_name)
  ),
  layout_columns(
    card(
      title = "Gene Structure",
      plotOutput("ggtranscript_plot")
    ),
    card(
      title = "Abundance",
      plotOutput("abundance")
    ),
    col_widths = c(8, 4)
  )
)

plot_ggtranscript <- function(gtf, gene_of_interest) {
  exons <- gtf %>%
    filter(gene_name == {{ gene_of_interest }})

  exons %>%
    ggplot(
      aes(xstart = start, xend = end, y = transcript_id, strand = strand)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id"),
    ) +
    geom_range(aes(fill = exon_status)) +
    facet_grid(dataset ~ ., scales = "free", space = "free") +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text.y = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "top"
    )
}

plot_SFARI <- function(gene_of_interest) {
  pbid_abundance_gene <- pbid_abundance %>%
    filter(transcript_id %in% (combined_gtf %>% filter(gene_name == {{ gene_of_interest }} & dataset == "SFARI") %>% pull("transcript_id"))) %>%
    pivot_longer(cols = c("CN", "NPC", "iPSC"), names_to = "time_point", values_to = "abundance")

  pbid_abundance_gene$time_point <- factor(pbid_abundance_gene$time_point, levels = c("iPSC", "NPC", "CN"))

  pbid_abundance_gene %>%
    ggplot(aes(x = transcript_id, y = abundance, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("SFARI") +
    ylab("log2(CPM + 1)") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    coord_flip() +
    scale_fill_manual(values = c(iPSC = "#e0a19c", NPC = "#e14bd5", CN = "#ff0011"))
}

plot_Patowary <- function(gene_of_interest) {
  talon_abundance %>%
    filter(transcript_id %in% (combined_gtf %>% filter(gene_name == {{ gene_of_interest }} & dataset == "Patowary et al.") %>% pull("transcript_id"))) %>%
    pivot_longer(cols = c("CP_mean", "VZ_mean"), names_to = "time_point", values_to = "abundance") %>%
    ggplot(aes(x = transcript_id, y = abundance, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("Patowary et al.") +
    ylab("log2(CPM + 1)") +
    scale_color_brewer(palette = "PuOr") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_manual(values = c(CP_mean = "#9797ea", VZ_mean = "#08a0a2")) +
    coord_flip()
}

plot_overlapped <- function(gene_of_interest) {
  transcript_id_to_plot <- combined_gtf %>%
    filter(gene_name == {{ gene_of_interest }} & (dataset == "Overlapped")) %>%
    pull(transcript_id) %>%
    unique()

  pbid_abundance_gene <- pbid_abundance %>%
    filter(transcript_id %in% transcript_id_to_plot) %>%
    rowwise() %>%
    summarise(
      transcript_id,
      mean_value = mean(c_across(CN:iPSC), na.rm = TRUE),
      dataset = "SFARI"
    ) %>%
    ungroup()

  talon_abundance_gene <- talon_abundance %>%
    filter(transcript_id %in% transcript_id_to_plot) %>%
    rowwise() %>%
    summarise(
      transcript_id,
      mean_value = mean(c_across(CP_mean:VZ_mean), na.rm = TRUE),
      dataset = "Patowary et al."
    ) %>%
    ungroup()

  bind_rows(pbid_abundance_gene, talon_abundance_gene) %>%
    ggplot(aes(x = transcript_id, y = mean_value, fill = dataset)) +
    geom_col(position = "dodge") +
    xlab("Overlapped") +
    ylab("log2(CPM + 1)") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_manual(values = c("SFARI" = "#f51c0c", "Patowary et al." = "#2908a2")) +
    coord_flip()
}

get_n_transcript_ratio <- function(gene_of_interest) {
  combined_gtf %>%
    filter(gene_name == {{ gene_of_interest }}) %>%
    group_by(dataset) %>%
    summarise(n_transcript_id = n_distinct(transcript_id)) %>%
    pull(n_transcript_id, name = dataset)
}

plot_abundance <- function(gene_of_interest) {
  if (setequal(c("Overlapped", "SFARI", "Patowary et al."), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_overlapped(gene_of_interest) / plot_Patowary(gene_of_interest) / plot_SFARI(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("SFARI"), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_SFARI(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("Patowary et al."), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_Patowary(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("Overlapped"), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_overlapped(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("SFARI", "Patowary et al."), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_Patowary(gene_of_interest) / plot_SFARI(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("SFARI", "Overlapped"), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_overlapped(gene_of_interest) / plot_SFARI(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  } else if (setequal(c("Patowary et al.", "Overlapped"), names(get_n_transcript_ratio(gene_of_interest)))) {
    plot_overlapped(gene_of_interest) / plot_Patowary(gene_of_interest) + plot_layout(heights = get_n_transcript_ratio(gene_of_interest) * 2)
  }
}

server <- function(input, output, session) {
  output$ggtranscript_plot <- renderPlot(
    {
      plot_ggtranscript(combined_gtf, input$select_genes)
    },
    height = reactive({
      20 * sum(get_n_transcript_ratio(input$select_genes)) + 100
    })
  )
  output$abundance <- renderPlot(
    {
      plot_abundance(input$select_genes)
    },
    height = reactive({
      20 * sum(get_n_transcript_ratio(input$select_genes)) + 100
    })
  )
}

shinyApp(ui, server)
