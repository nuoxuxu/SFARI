library(ggplot2)
library(ggtranscript)
library(dplyr)
library(readr)
library(here)

combined_gtf <- read_csv(here("proc/shiny.csv"))

plot_ggtranscript <- function(gtf, gene_of_interest) {
  exons <- gtf %>% 
    filter(feature == "exon") %>% 
    filter(gene_name == gene_of_interest)
  
  exons %>%
    ggplot(
      aes(xstart = start, xend = end, y = transcript_name, strand = strand)
    ) +
    geom_range() +
    geom_intron(
      data = to_intron(exons, "transcript_name"),
    ) +
    facet_grid(dataset ~ ., scales = "free", space = "free") +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
    )
}

server <- function(input, output, session) {
  updateSelectizeInput(session, "select_genes", choices = unique(combined_gtf$gene_name), server = TRUE)
  output$ggtranscript_plot <- renderPlot({
    plot_ggtranscript(combined_gtf, input$select_genes)
  })
}