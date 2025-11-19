library(ggplot2)
library(ggtranscript)
library(dplyr)
library(shiny)
library(bslib)
library(patchwork)
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

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00",
    "Other" = "#000000"
)

pastelColorVector <- c(
    "FSM" = "#7FCEB9",
    "ISM" = "#7FB8D8",
    "NIC" = "#EAAE7F",
    "NNC" = "#F2CF7F",
    "Other" = "#7F7F7F"
)

gtf <- read.csv("data/gtf.csv")
lr_log2_cpm <- read.csv("data/lr_log2_cpm.csv") %>% 
  mutate(
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC")),
    time_point = factor(time_point, levels = c("t30", "t04", "t00"))
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
      data = to_intron(exons, "transcript_id")
    ) +
    # hide exon legend
    geom_range(
      aes(fill = structural_category),
      height = 0.25,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = pastelColorVector, guide = "none") +
    ggnewscale::new_scale_fill() +
    # keep only CDS legend
    geom_range(
      data = CDS,
      aes(fill = structural_category),
      color = "black"
    ) +
    scale_fill_manual("Structural Category", values = colorVector) +
    facet_grid(structural_category ~ ., scales = "free", space = "free")
}

plot_abundance <- function(gene_of_interest) {
  lr_log2_cpm %>% 
    filter(associated_gene == gene_of_interest) %>%
    ggplot(aes(x = abundance, y = isoform, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("mean log2(CPM + 1)") +
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
  gtf %>%
    filter(type == "transcript") %>%
    filter(associated_gene == gene_of_interest) %>%
    nrow()
}

ui <- fluidPage(
  titlePanel("Human iPSC-derived Neurons Full-length Transcriptome"),
  helpText("Explore transcript structures and abundances for genes of interest. Transcript structural categories are based on comparison to GENCODE v47."),
  sidebarPanel(
    selectizeInput(
      "select_genes",
      "Select gene of interest:",
      choices = choices <- unique(gtf$associated_gene),
      selected = "AGO1"
    ),
    helpText("To view the genomic region of the selected gene in the UCSC Genome Browser, click the link below:"),
    #insert a hyperlink here
    uiOutput("ucsc_link"),
    width = 3
  ),
  mainPanel(
    plotOutput("combined_plot")
  )
)

server <- function(input, output, session) {
  output$combined_plot <- renderPlot(
    {
      validate(need(input$select_genes, 'Choose a gene!'))
      plot_combined(input$select_genes)
    },
    height = reactive({
      35 * sum(get_n_transcripts(input$select_genes)) + 100
    })
  )
  output$ucsc_link <- renderUI({
    validate(need(input$select_genes, 'Choose a gene!'))
  gene_info <- gtf %>%
    filter(type == "transcript", associated_gene == input$select_genes)
  chrom <- gene_info %>% select(seqnames) %>% distinct() %>% pull()
  start <- gene_info %>% select(start) %>% min()
  end <- gene_info %>% select(end) %>% max()
    url <- paste0("https://genome.ucsc.edu/s/nuoxuxu/iPSC_Neuron_Proteogenomic_atlas?position=", chrom, ":", start, "-", end)
    tags$a(href = url, "View in UCSC Genome Browser", target = "_blank")
  })
}

shinyApp(ui, server)