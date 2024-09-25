library(ggplot2)
library(ggtranscript)
library(dplyr)
library(readr)
library(here)

combined_gtf <- read_csv(here("./shiny.csv"))

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

server <- function(input, output, session) {
  updateSelectizeInput(session, "select_genes", choices = unique(combined_gtf$gene_name), server = TRUE)
  output$ggtranscript_plot <- renderPlot(
    {plot_ggtranscript(combined_gtf, input$select_genes)},
    height = reactive({20 * get_n_transcripts(combined_gtf, input$select_genes) + 100})
    )
    }