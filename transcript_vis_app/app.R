library(ggplot2)
library(ggtranscript)
library(dplyr)
library(tidyr)
library(readr)
library(bslib)

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
    layout_columns(
      card(
        card_header("Overlapped transcripts abundance"),
        uiOutput("Overlapped_abundance"),
      ),
      card(
        card_header("Patowary et al. transcripts abundance"),
        uiOutput("Patowary_abundance"),
      ),
      card(
        card_header("SFARI transcripts abundance"),
        uiOutput("SFARI_abundance")
      ),
      col_widths = c(12, 12, 12)
    ),
    col_widths = c(8, 4)
  )
)

plot_ggtranscript <- function(gtf, gene_of_interest) {
  exons <- gtf %>% 
    filter(gene_name == {{gene_of_interest}})
  
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
        axis.text.y = element_text(size = 12)
    )
}

get_n_transcripts <- function(gtf, gene_of_interest) {
  gtf %>%
    filter(gene_name == {{gene_of_interest}}) %>% 
    summarise(n_distinct(transcript_id)) %>% 
    pull()
    }

plot_SFARI <- function(gene_of_interest) {
  pbid_abundance_gene <- pbid_abundance %>% 
    filter(transcript_id %in% (combined_gtf %>% filter(gene_name == {{gene_of_interest}}) %>% pull("transcript_id"))) %>% 
    pivot_longer(cols = c("CN", "NPC", "iPSC"), names_to = "time_point", values_to = "abundance")

  pbid_abundance_gene$time_point <- factor(pbid_abundance_gene$time_point, levels = c("iPSC", "NPC", "CN"))  
    
  gg <- pbid_abundance_gene %>% 
    ggplot(aes(x=transcript_id, y=abundance, fill=time_point)) +
    geom_col(position="dodge") +
    ylab("counts per million (CPM)") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank()
    ) + coord_flip()
  renderPlot(gg)
}

plot_Patowary <- function(gene_of_interest) {
  gg <- talon_abundance %>% 
    filter(transcript_id %in% (combined_gtf %>% filter(gene_name == {{gene_of_interest}}) %>% pull("transcript_id"))) %>% 
    pivot_longer(cols = c("CP_mean", "VZ_mean"), names_to = "time_point", values_to = "abundance") %>%
    ggplot(aes(x=transcript_id, y=abundance, fill=time_point)) +
    geom_col(position="dodge") +
    ylab("counts per million (CPM)") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x = element_blank()
    ) + coord_flip()
  renderPlot(gg)
}

plot_overlapped <- function(gene_of_interest) {
  transcript_id_to_plot <- combined_gtf %>%
    filter(gene_name == {{gene_of_interest}} & (dataset == "Overlapped")) %>%
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

  gg <- bind_rows(pbid_abundance_gene, talon_abundance_gene) %>% 
  ggplot(aes(x=transcript_id, y=mean_value, fill=dataset)) +
  geom_col(position="dodge") +
  ylab("counts per million (CPM)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank()
  ) + coord_flip()
  renderPlot(gg)
}

Overlapped_no_transcripts <- function(gene_of_interest) {
  combined_gtf %>%
    filter(gene_name == {{gene_of_interest}} & dataset == "Overlapped") %>%
    pull(transcript_id) %>%
    length()
}

TALON_no_transcripts <- function(gene_of_interest) {
  combined_gtf %>%
    filter(gene_name == {{gene_of_interest}} & dataset == "Patowary et al.") %>%
    pull(transcript_id) %>%
    length()
}

SFARI_no_transcripts <- function(gene_of_interest) {
  combined_gtf %>%
    filter(gene_name == {{gene_of_interest}} & dataset == "SFARI") %>%
    pull(transcript_id) %>%
    length()
}

server <- function(input, output, session) {
  Overlapped_show <- reactive({Overlapped_no_transcripts(input$select_genes) != 0})
  TALON_show <- reactive({TALON_no_transcripts(input$select_genes) != 0})
  SFARI_show <- reactive({SFARI_no_transcripts(input$select_genes) != 0})

  output$ggtranscript_plot <- renderPlot(
    {
      plot_ggtranscript(combined_gtf, input$select_genes)
      },
    height = reactive({20 * get_n_transcripts(combined_gtf, input$select_genes) + 100})
  )
  output$Overlapped_abundance <- renderUI(
    if (Overlapped_show()) {
      plot_overlapped(input$select_genes)
    }
  )
  output$Patowary_abundance <- renderUI(
    if (TALON_show()) {
      plot_Patowary(input$select_genes)
    }
  )  
  output$SFARI_abundance <- renderUI(
    if (SFARI_show()) {
      plot_SFARI(input$select_genes)
    }
  )
}

shinyApp(ui, server)