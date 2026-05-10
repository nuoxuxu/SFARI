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

time_point_colors <- c("t00" = "#e0a19c", "t04" = "#e14bd5", "t30" = "#ff0011")

gtf <- read.csv("data/gtf.csv")
lr_log2_cpm <- read.csv("data/lr_log2_cpm.csv") %>%
  mutate(
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC")),
    time_point = factor(time_point, levels = c("t30", "t04", "t00")),
    has_CDS = as.logical(has_CDS)
  )

apply_coding_filter <- function(data, coding_filter, id_col = "isoform") {
  if (coding_filter == "coding") {
    data %>% filter(has_CDS == TRUE)
  } else if (coding_filter == "non_coding") {
    data %>% filter(has_CDS == FALSE)
  } else {
    data
  }
}

plot_ggtranscript <- function(gene_of_interest, coding_filter = "all") {
  gtf_gene <- gtf %>%
    filter(associated_gene == gene_of_interest) %>%
    apply_coding_filter(coding_filter)

  exons <- gtf_gene %>% filter(type == "exon")
  CDS <- gtf_gene %>% filter(type == "CDS")

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

plot_abundance <- function(gene_of_interest, coding_filter = "all") {
  data <- lr_log2_cpm %>%
    filter(associated_gene == gene_of_interest) %>%
    apply_coding_filter(coding_filter)

  data %>%
    ggplot(aes(x = log2_cpm, y = isoform, color = time_point)) +
    # Individual replicate points
    geom_point(
      position = position_dodge(width = 0.7),
      alpha = 0.5,
      size = 1.5
    ) +
    # Mean ± SEM
    stat_summary(
      fun.data = mean_se,
      geom = "pointrange",
      position = position_dodge(width = 0.7),
      size = 0.6,
      linewidth = 0.8
    ) +
    xlab("log2(CPM + 1)") +
    theme(
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 14),
    ) +
    scale_color_manual("Time point", breaks = c("t00", "t04", "t30"), values = time_point_colors) +
    facet_grid(structural_category ~ ., scales = "free", space = "free")
}

plot_combined <- function(gene_of_interest, coding_filter = "all") {
  ggtranscript_plot <- plot_ggtranscript(gene_of_interest, coding_filter)
  abundance <- plot_abundance(gene_of_interest, coding_filter)
  ggtranscript_plot | abundance + plot_layout(widths = c(2, 1))
}

get_n_transcripts <- function(gene_of_interest, coding_filter = "all") {
  gtf %>%
    filter(type == "transcript", associated_gene == gene_of_interest) %>%
    apply_coding_filter(coding_filter) %>%
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
    selectInput(
      "coding_filter",
      "Filter isoforms by coding status:",
      choices = c(
        "All isoforms"              = "all",
        "Protein-coding (has CDS)"  = "coding",
        "Non-coding (no CDS)"       = "non_coding"
      ),
      selected = "all"
    ),
    helpText("To view the genomic region of the selected gene in the UCSC Genome Browser, click the link below:"),
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
      plot_combined(input$select_genes, input$coding_filter)
    },
    height = reactive({
      35 * sum(get_n_transcripts(input$select_genes, input$coding_filter)) + 100
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