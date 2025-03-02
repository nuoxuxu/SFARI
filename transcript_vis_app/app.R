library(ggplot2)
library(ggtranscript)
library(dplyr)
library(tidyr)
library(readr)
library(shiny)
library(bslib)
library(patchwork)

combined_gtf <- read_csv("./shiny.csv")
pbid_abundance <- read_csv("./pbid_abundance.csv")
talon_abundance <- read_csv("./talon_abundance.csv")
CDS_gtf <- read_csv("genome_gff3_gtf.csv")
gencode_gtf <- read_csv("GENCODE_v39.csv")
peptides_gtf <- read_csv("SFARI_peptides.csv")

ui <- page_fillable(
  navset_tab(
    nav_panel(
      title = "Comparison",
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
    ),
    nav_panel(
      title = "CDS",
        selectizeInput(
          "select_genes_CDS",
          "Select genes from list below:",
          choices = choices <- unique(CDS_gtf$gene_name)
        ),
        selectizeInput(
          "select_protein_class",
          "Select SQANTI protrein classes from list below:",
          choices = choices <- c(unique(as.character(unique(CDS_gtf[["protein_classification_base"]]))), "all"),
          selected = "all"
        ),        
        selectizeInput(
          "select_mapped_to_nov_trans_only",
          "Filter for peptides that represent novel splice junctions or mono-exons",
          choices = choices <- c("known", "novel", "all"),
          selected = "all"
      ),
      card(
        title = "Unique CDS",
        plotOutput("CDS_plot")
      )
    )
  )
)

plot_ggtranscript <- function(gene_of_interest) {
  exons <- combined_gtf %>%
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
    pivot_longer(cols = c("iPSC", "NPC", "CN"), names_to = "time_point", values_to = "abundance")

  pbid_abundance_gene$time_point <- factor(pbid_abundance_gene$time_point, levels = c("CN", "NPC", "iPSC"))

  pbid_abundance_gene %>%
    ggplot(aes(x = abundance, y = transcript_id, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("log2(CPM + 1)") +
    ylab("SFARI") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_manual(breaks = c("iPSC", "NPC", "CN"), values = c("#e0a19c", "#e14bd5", "#ff0011"))
}

plot_Patowary <- function(gene_of_interest) {
  talon_abundance %>%
    filter(transcript_id %in% (combined_gtf %>% filter(gene_name == {{ gene_of_interest }} & dataset == "Patowary et al.") %>% pull("transcript_id"))) %>%
    pivot_longer(cols = c("CP_mean", "VZ_mean"), names_to = "time_point", values_to = "abundance") %>%
    ggplot(aes(x = abundance, y = transcript_id, fill = time_point)) +
    geom_col(position = "dodge") +
    xlab("log2(CPM + 1)") +
    ylab("Patowary et al.") +
    scale_color_brewer(palette = "PuOr") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_manual(breaks = c("CP_mean", "VZ_mean"), values = c("#9797ea", "#08a0a2"))
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
    ggplot(aes(x = mean_value, y = transcript_id, fill = dataset)) +
    geom_col(position = "dodge") +
    xlab("log2(CPM + 1)") +
    ylab("Overlapped") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_fill_manual(breaks = c("SFARI", "Patowary et al."), values = c("#f51c0c", "#2908a2"))
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

plot_gencode <- function(gtf, gene_of_interest, xmin, xmax) {
    exons <- gtf %>%
        filter(gene_name == {{ gene_of_interest }}) %>%
        filter(feature %in% c("exon"))
    CDS <- gtf %>%
        filter(gene_name == {{ gene_of_interest }}) %>%
        filter(feature == "CDS")

    exons %>%
        ggplot(
            aes(xstart = start, xend = end, y = transcript_id, strand = strand)
            ) +
        geom_range(
            fill = "white",
            height = 0.25
            ) +            
        geom_intron(
            data = to_intron(exons, "transcript_id")
            ) +
        geom_range(
            data = CDS
            ) +
        coord_cartesian(xlim = c({{xmin}}, {{xmax}})) +
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

plot_CDS <- function(gtf, gene_of_interest, select_protein_class) {
    if ((select_protein_class == "all")) {
      gtf <- gtf %>%
          filter(gene_name == {{gene_of_interest}})
    } else {
      gtf <- gtf %>%
          filter(gene_name == {{gene_of_interest}}) %>%
          filter(protein_classification_base == {{select_protein_class}})
    }

    gtf %>%
        ggplot(
            aes(xstart = start, xend = end, y = transcript_id, strand = strand)
            ) +
        geom_range(
            aes(fill = protein_classification_base)
            ) +  
        geom_intron(
            data = to_intron(gtf, "transcript_id")
            ) +
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

plot_peptides <- function(gtf, gene_of_interest, novelty, xmin, xmax) {
    if (novelty == "all") {
      gtf <- gtf %>%
          filter(gene_name == {{gene_of_interest}})
    } else {
      gtf <- gtf %>%
          filter(gene_name == {{gene_of_interest}})%>%
          filter(novelty == {{novelty}})
    }
    gtf %>%
        ggplot(
            aes(xstart = start, xend = end, y = transcript_id, strand = strand)
            ) +
        geom_range(
            aes(fill = novelty)
            ) +  
        geom_intron(
            data = to_intron(gtf, "transcript_id")
            ) +
        coord_cartesian(xlim = c({{xmin}}, {{xmax}})) +
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

get_n_transcript_CDS <- function(gene_of_interest, select_protein_class, novelty) {
    if (select_protein_class == "all") {
      n_CDS <- CDS_gtf %>%
          filter(gene_name == {{gene_of_interest}}) %>%
          summarise(CDS = n_distinct(transcript_id))
    } else {
      n_CDS <- CDS_gtf %>%
          filter(gene_name == {{gene_of_interest}}) %>%
          filter(protein_classification_base == {{select_protein_class}}) %>%
          summarise(CDS = n_distinct(transcript_id))
    }
    if (novelty == "all") {
      n_peptides <- peptides_gtf %>%
          filter(gene_name == {{gene_of_interest}}) %>%
          summarise(peptides = n_distinct(transcript_id))
    } else {
      n_peptides <- peptides_gtf %>%
          filter(gene_name == {{gene_of_interest}}) %>%
          filter(novelty == {{novelty}}) %>%
          summarise(peptides = n_distinct(transcript_id))      
    }
    n_gencode <- gencode_gtf %>%
        filter(gene_name == {{gene_of_interest}}) %>%
        summarise(gencode = n_distinct(transcript_id))    

    cbind(n_gencode, n_CDS, n_peptides)
}

plot_combined <- function(gene_of_interest, select_protein_class, select_mapped_to_nov_trans_only) {
    CDS <- plot_CDS(CDS_gtf, {{gene_of_interest}}, {{select_protein_class}})
    xmin <- ggplot_build(CDS)$layout$panel_params[[1]]$x.range[1]
    xmax <- ggplot_build(CDS)$layout$panel_params[[1]]$x.range[2]
    gencode <- plot_gencode(gencode_gtf, {{gene_of_interest}}, {{xmin}}, {{xmax}})
    peptides <- plot_peptides(peptides_gtf, {{gene_of_interest}}, {{select_mapped_to_nov_trans_only}}, {{xmin}}, {{xmax}})
    gencode / CDS / peptides + plot_layout(heights = get_n_transcript_CDS(gene_of_interest, select_protein_class, select_mapped_to_nov_trans_only) * 2)
}

server <- function(input, output, session) {
  output$ggtranscript_plot <- renderPlot(
    {
      plot_ggtranscript(input$select_genes)
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
  output$CDS_plot <- renderPlot(
    {
      plot_combined(input$select_genes_CDS, input$select_protein_class, input$select_mapped_to_nov_trans_only)
    },
    height = reactive({
      20 * sum(get_n_transcript_CDS(input$select_genes_CDS, input$select_protein_class, input$select_mapped_to_nov_trans_only)) + 100
    })
  )
}

shinyApp(ui, server)
