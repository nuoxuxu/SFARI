### MAKING FIGURES ###
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(readxl)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(arrow)
library(gprofiler2)



setwd("")


# FIGURE 3A ---------------------------------------------------------------
sample_data <- data.frame(timepoint=c("t00", "t04", "t30", 
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30",
                                      "t00", "t04", "t30", 
                                      "t00", "t04", "t30"),
                          expr=c(2.5, 2.5, 2.5, 
                                 2.5, 2.5, 0.5, 
                                 2.5, 2.5, 4.5, 
                                 4.5, 2.5, 2.5, 
                                 4.5, 2.5, 0.5, 
                                 4, 2, 4, 
                                 0.5, 2.5, 2.5, 
                                 1, 3, 1, 
                                 0.5, 2.5, 4.5), 
                          cluster = c("--", "--", "--",
                                      "-D", "-D", "-D",
                                      "-U", "-U", "-U", 
                                      "D-", "D-", "D-", 
                                      "DD", "DD", "DD", 
                                      "DU", "DU", "DU",
                                      "U-", "U-", "U-",
                                      "UD", "UD", "UD",
                                      "UU", "UU", "UU"))
sample_data$timepoint <- as.factor(sample_data$timepoint)
sample_data$cluster <- as.factor(sample_data$cluster)



safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


Fig_3A <- ggplot(data=sample_data, aes(x=timepoint, y=expr, group = 1)) +
  geom_line(aes(colour = cluster), linewidth = 1.5) +
  geom_point(aes(colour = cluster), size = 2.75) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5)) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.position = "none", 
        strip.background = element_blank(), 
        strip.text = element_text(face="bold", size=8)) +
  scale_colour_manual(values = safe_colorblind_palette[1:9]) +
  facet_wrap(sample_data$cluster, 3, 3, axes = "all") +
  labs(x = "Time", y = "RNA or Protein\nExpression")
Fig_3A



# RIVER PLOT --------------------------------------------------------------
# GENES VS PROTEINS
#Compare genes and proteins. Riverplots each cluster vs each cluster.
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

gene_clusters <- read.csv("./code/expression/output/gene_clusters_v2.csv")
protein_clusters <- read.csv("./code/expression/output/protein_clusters_v3.csv")


genesum_prot_clusters <- gene_clusters[,c(1,11)]
colnames(genesum_prot_clusters) <- c("gene_name", "mRNA")

genesum_prot_clusters <- merge(genesum_prot_clusters, protein_clusters[, c("GeneSymbol", "cluster")], 
                               by.x = "gene_name", by.y = "GeneSymbol", all.x = TRUE)

genesum_prot_clusters <- na.omit(genesum_prot_clusters) #8493 remaining out of 8878
colnames(genesum_prot_clusters)[3] <- "Protein"

genesum_prot_long <- genesum_prot_clusters %>% pivot_longer(cols = c(mRNA, Protein),
                                                            names_to = "gene_prot",
                                                            values_to = "cluster")


genesum_prot_long$cluster <- as.factor(genesum_prot_long$cluster)
genesum_prot_long$rep <- rep(1:nrow(genesum_prot_clusters), each=2)

label_fix <- data.frame(
  gene_prot = c("Protein", "Protein"),
  y = c(1131.5, 1012.5),            
  label = c("DD", "DU"), 
  x = c(1.95, 2.05))

genesum_prot_long$label <- as.character(genesum_prot_long$cluster)
genesum_prot_long$label[genesum_prot_long$cluster %in% c("DD", "DU") & 
                          genesum_prot_long$gene_prot == "Protein"] <- ""


genesum_river <- ggplot(genesum_prot_long,
                        aes(x = gene_prot, stratum = cluster, alluvium = rep,
                            fill = cluster, label = label)) +
  geom_alluvium(aes(fill = cluster),
                curve_type = "sigmoid") +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = label), 
            fontface = "bold", 
            size = 3.75) +
  scale_fill_manual(values = safe_colorblind_palette[1:9]) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8600)) +
  scale_x_discrete(expand = c(0, 0)) + #move y axis closer
  theme_void(base_size = 8) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text = element_text(size = 8),
        axis.ticks.y = element_line(size = 1),
        axis.ticks.length.y = unit(3, "pt"),
        axis.line.y = element_line(size = 0.75, colour = "black"), 
        plot.margin = margin(
          t = 6,
          r = 6,
          b = 6,
          l = 6,
          unit = "pt" ))

genesum_river_updated <- genesum_river +
  geom_text(data = label_fix,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            fontface = "bold", size = 3.75)
genesum_river_updated


#
table(genesum_prot_clusters$mRNA) #59% changing
table(genesum_prot_clusters$Protein) #42% changing



# EXAMPLE GENES -----------------------------------------------------------
###Assigning concordant (0) and discordant 1 or 2
genesum_prot_clusters$mRNA <- as.character(genesum_prot_clusters$mRNA)
genesum_prot_clusters$Protein <- as.character(genesum_prot_clusters$Protein)

genesum_prot_clusters[,c("Genes_1", "Genes_2")] <- str_split_fixed(genesum_prot_clusters$mRNA, "", 2)
genesum_prot_clusters[,c("Proteins_1", "Proteins_2")] <- str_split_fixed(genesum_prot_clusters$Protein, "", 2)

genesum_prot_clusters$match1 <- ifelse((genesum_prot_clusters$Genes_1 == genesum_prot_clusters$Proteins_1), 0, 1)
genesum_prot_clusters$match2 <- ifelse((genesum_prot_clusters$Genes_2 == genesum_prot_clusters$Proteins_2), 0, 1)
genesum_prot_clusters$Level <- genesum_prot_clusters$match1 + genesum_prot_clusters$match2

sum(genesum_prot_clusters$Level == 0)/nrow(genesum_prot_clusters) #45.5
sum(genesum_prot_clusters$Level == 1)/nrow(genesum_prot_clusters) #42.5
sum(genesum_prot_clusters$Level == 2)/nrow(genesum_prot_clusters) #12.0


#Intersect the level2 with SFARI genes
SFARI <- read.csv("./data/SFARI-Gene_genes_04-03-2025release_04-15-2025export.csv")

level_SFARI <- genesum_prot_clusters[genesum_prot_clusters$gene_name %in% SFARI$gene.symbol, ]
level_SFARI <- level_SFARI %>% filter((mRNA != "--" | Protein != "--"))


#Make example plots for gene to protein fig. 
#Prepare gene data. 
gene_cpm <- read.csv("./code/expression/output/gene_CPM.csv")

gene_SFARI <- gene_cpm[gene_cpm$gene_name %in% level_SFARI$gene_name,]
gene_SFARI <- merge(gene_SFARI, gene_clusters[, c("gene_name", "cluster")],
                     by.x = "gene_name", by.y = "gene_name", all.x = TRUE)

#Pivot longer.
gene_summ <- gene_SFARI %>% pivot_longer(
  cols = starts_with("sum_"), # Select columns that start with "sum_"
  names_to = c("Type", "Group"), # Create two new columns: "Type" (iPSC, NPC, CN) and "Group"
  names_sep = "_", # Split the names at "_"
  values_to = "Value")

#Log2 transform the data, then calculate mean and 95% CI, and SD
gene_long <- gene_summ %>%
  mutate(Value_log2 = log2(Value + 1))

gene_long$Group <- factor(gene_long$Group, levels = c("iPSC", "NPC", "CN"),
                        labels = c("t00", "t04", "t30"))

# Calculate mean and 95% CI, and SD
gene_summary <- gene_long %>%
  group_by(gene_name, Group) %>%
  summarise(
    Mean = mean(Value_log2),
    SE = sd(Value_log2) / sqrt(n()),
    CI_Lower = Mean - qt(0.975, df = n() - 1) * SE,
    CI_Upper = Mean + qt(0.975, df = n() - 1) * SE, 
    sd = sd(Value_log2, na.rm = TRUE) )
gene_summary <- merge(gene_summary, gene_clusters[, c("gene_name", "cluster")],
                    by.x = "gene_name", by.y = "gene_name", all.x = F)

###Prepare protein data
prot_quant <- read.csv("./code/expression/output/MSstatsTMT_abundance.csv")
final_filtered_diff_prot <- read.csv("./code/expression/output/MSstatsTMT_diff_prot_expr_filtered.csv")

#Calculate median across fractions:
prot_quant_median <- prot_quant %>% group_by(Protein, Channel) %>%
  summarise(Median_Abundance = median(Abundance)) %>%
  left_join(prot_quant %>% dplyr::select(Channel, Condition) %>% distinct(Channel, Condition) , by = "Channel")

summary_prot_abund <- prot_quant_median %>%
  group_by(Protein, Condition) %>%
  summarise(
    mean = mean(Median_Abundance, na.rm = TRUE),
    sd = sd(Median_Abundance, na.rm = TRUE),
    n = n()
  )

summary_prot_abund <- summary_prot_abund %>% 
  left_join(final_filtered_diff_prot %>% dplyr::select(Protein, GeneSymbol) %>% distinct(Protein, GeneSymbol), by = "Protein")
summary_prot_abund <- summary_prot_abund %>% left_join(genesum_prot_clusters %>% dplyr::select(gene_name, Protein), by = c("GeneSymbol" = "gene_name"))

summary_prot_abund$Condition <- factor(summary_prot_abund$Condition, levels = c("t00", "t4", "t30"), 
                                       labels = c("t00", "t04", "t30"))


#Chosen examples
gene_ex <- gene_summary[gene_summary$gene_name %in% level_SFARI$gene_name,]
prot_ex <- summary_prot_abund[summary_prot_abund$GeneSymbol %in% level_SFARI$gene_name,]

#Make 6 indiv plots per gene b/c I need gene name on the left, and specific scales
THBS1_gene <- gene_ex %>% filter(gene_name == "THBS1") %>% ggplot(aes(x = Group, y = Mean, group = gene_name, color = cluster)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(0, 4.2)) +
  theme_classic(base_size = 8) +
  labs(title = "mRNA", subtitle = "log2(CPM + 1)", y = "THBS1") +
  scale_color_manual(values=safe_colorblind_palette[4], 
                     breaks = c("D-")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5, size = 6.4), 
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold")) +
  geom_text(aes(x=Inf,y=Inf,
                hjust=c(1.5),
                vjust=c(1.5),
                label="D-",
                fontface = "bold"), 
            size = 4)
THBS1_gene

MYO1E_gene <- gene_ex %>% filter(gene_name == "MYO1E") %>% ggplot(aes(x = Group, y = Mean, group = gene_name, color = cluster)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(1, 8)) +
  theme_classic(base_size = 8) +
  labs(y = "MYO1E") +
  scale_color_manual(values=safe_colorblind_palette[4], 
                     breaks = c("D-")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold")) +
  geom_text(aes(x=Inf,y=Inf,hjust=c(1.5),
                vjust=c(1.5),label="D-", fontface = "bold"), 
            size = 4)
MYO1E_gene

DPP4_gene <- gene_ex %>% filter(gene_name == "DPP4") %>% ggplot(aes(x = Group, y = Mean, group = gene_name, color = cluster)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(-0.1, 4)) +
  theme_classic(base_size = 8) +
  labs(y = "DPP4") +
  scale_color_manual(values=safe_colorblind_palette[4], 
                     breaks = c("D-")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold")) +
  geom_text(aes(x=Inf,y=Inf,hjust=c(1.5),
                vjust=c(1.5),label="D-", fontface = "bold"), 
            size = 4)
DPP4_gene

#Proteins
THBS1_prot <- prot_ex %>% filter(GeneSymbol == "THBS1") %>% ggplot(aes(x = Condition, y = mean, group = Protein.y, color = Protein.y)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(3.8, 7.5)) +
  theme_classic(base_size = 8) +
  labs(title = "Protein", subtitle = "log2(Abundance)") +
  scale_color_manual(values=safe_colorblind_palette[c(4)], 
                     breaks = c("D-")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), 
        plot.subtitle = element_text(hjust = 0.5, size = 6.4), 
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title = element_blank()) +
  geom_text(aes(x=Inf,y=Inf,hjust=c(1.5),
                vjust=c(1.5),label="D-", fontface = "bold"), 
            size = 4)
THBS1_prot

MYO1E_prot <- prot_ex %>% filter(GeneSymbol == "MYO1E") %>% ggplot(aes(x = Condition, y = mean, group = Protein.y, color = Protein.y)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(3.8, 7.5)) +
  theme_classic(base_size = 8) +
  scale_color_manual(values=safe_colorblind_palette[5], 
                     breaks = c("DD")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title = element_blank()) +
  geom_text(aes(x=Inf,y=Inf,hjust=c(1.5),
                vjust=c(1.5),label="DD", fontface = "bold"), 
            size = 4)
MYO1E_prot

DPP4_prot <- prot_ex %>% filter(GeneSymbol == "DPP4") %>% ggplot(aes(x = Condition, y = mean, group = Protein.y, color = Protein.y)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2.75) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, linewidth = 1) +
  scale_y_continuous(limits = c(3.8, 7.5)) +
  theme_classic(base_size = 8) +
  scale_color_manual(values=safe_colorblind_palette[c(4,1,8)], 
                     breaks = c("D-", "--", "UD")) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.position = "none", 
        axis.text = element_text(size = 8),
        axis.title = element_blank()) +
  geom_text(aes(x=Inf,y=Inf,hjust=c(1.5),
                vjust=c(1.5),label="UD", fontface = "bold"), 
            size = 4)
DPP4_prot


combined_plot <- ggarrange(THBS1_gene, THBS1_prot, 
                           MYO1E_gene, MYO1E_prot, 
                           DPP4_gene, DPP4_prot ,
                           ncol = 2,
                           nrow = 3, 
                           heights = c(1.32, 1, 1))
combined_plot


# 9 CLUSTER ASD RISK GENE ENRICHMENT------------------------------------------
####GO term enrichment for 9 protein clusters:
#bg proteins is all genes with mRNA and protein expression
bg.prot <- genesum_prot_clusters %>% pull(gene_name) %>% unique()


#use 0 to replace -
gene_lists <- list(genes_00 = genesum_prot_clusters %>% filter(Protein == "--") %>% pull(gene_name) %>% unique(),
                   genes_0D = genesum_prot_clusters %>% filter(Protein == "-D") %>% pull(gene_name) %>% unique(),
                   genes_U0 = genesum_prot_clusters %>% filter(Protein == "U-") %>% pull(gene_name) %>% unique(),
                   genes_0U = genesum_prot_clusters %>% filter(Protein == "-U") %>% pull(gene_name) %>% unique(),
                   genes_D0 = genesum_prot_clusters %>% filter(Protein == "D-") %>% pull(gene_name) %>% unique(),
                   genes_UU = genesum_prot_clusters %>% filter(Protein == "UU") %>% pull(gene_name) %>% unique(),
                   genes_UD = genesum_prot_clusters %>% filter(Protein == "UD") %>% pull(gene_name) %>% unique(),
                   genes_DU = genesum_prot_clusters %>% filter(Protein == "DU") %>% pull(gene_name) %>% unique(),
                   genes_DD = genesum_prot_clusters %>% filter(Protein == "DD") %>% pull(gene_name) %>% unique())

###Run GO term search:
run_GO_enrichment <- function(gene_list, bg, searches) {
  go_result <- gost(gene_list, 
                    organism = "hsapiens", 
                    significant = T,
                    exclude_iea = F, 
                    user_threshold = 0.05,
                    correction_method = "fdr", 
                    custom_bg = bg, 
                    sources = searches)
  return(go_result)
}
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = bg.prot,
                                                                  searches = searches)
}

enrichment_results$genes_00$result$query <- "--"
enrichment_results$genes_0D$result$query <- "-D"
enrichment_results$genes_U0$result$query <- "U-"
enrichment_results$genes_0U$result$query <- "-U"
enrichment_results$genes_D0$result$query <- "D-"
enrichment_results$genes_UU$result$query <- "UU"
enrichment_results$genes_UD$result$query <- "UD"
enrichment_results$genes_DU$result$query <- "DU"
enrichment_results$genes_DD$result$query <- "DD"


GO_list <- rbind(enrichment_results$genes_00$result, 
                 enrichment_results$genes_0D$result, 
                 enrichment_results$genes_U0$result, 
                 enrichment_results$genes_0U$result, 
                 enrichment_results$genes_D0$result, 
                 enrichment_results$genes_UU$result, 
                 enrichment_results$genes_UD$result, 
                 enrichment_results$genes_DU$result, 
                 enrichment_results$genes_DD$result)

GO_list$query <- as.character(GO_list$query)
write.csv(GO_list[,-14], "./tables/protein_clusters_GO_terms.csv", 
          row.names = F)


#########SFARI Gene Enrichment
#SFARI genes subset to only proteins expressed:
SFARI_prot <- SFARI %>% filter(gene.symbol %in% genesum_prot_clusters$gene_name) #805

# Function to compute Fisher's exact test for pairwise comparisons
compute_fisher_enrich <- function(setA, setB, bg) {
  # Calculate overlap and unique genes
  overlap <- length(intersect(setA, setB))
  inAnotB <- length(setdiff(setA, setB))
  inBnotA <- length(setdiff(setB, setA))
  inNeither <- length(bg) - length(union(setA, setB))
  
  # Create the contingency table
  contingency_table <- matrix(c(overlap, inAnotB, inBnotA, inNeither), nrow = 2) 
  
  # Perform Fisher's Exact Test
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  return(c(overlap = overlap, inAnotB = inAnotB, inBnotA = inBnotA, inNeither = inNeither, p_value = fisher_result$p.value))
}


#Testing for enrichment in the each cluster:
results_list <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI_prot$gene.symbol
  setB <- gene_lists[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.prot)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists[i]))
  results_list[[pair_name]] <- result
}

results_df_SFARI_prot <- do.call(rbind, results_list)
results_df_SFARI_prot <- data.frame(pair = names(results_list), results_df_SFARI_prot)
results_df_SFARI_prot$p_adjusted <- p.adjust(results_df_SFARI_prot$p_value, method = "bonferroni")
results_df_SFARI_prot$gene <- "Protein"

SFARI_gene_enrich <- results_df_SFARI_prot

SFARI_gene_enrich$cluster <- gsub("0", "-", substr(SFARI_gene_enrich$pair, nchar(SFARI_gene_enrich$pair)-1, nchar(SFARI_gene_enrich$pair)))
SFARI_gene_enrich$cluser <- as.factor(SFARI_gene_enrich$cluster)
Fig3d_title <- "SFARI Gene Enrichment"

SFARI_gene_enrich$gene <- as.factor(SFARI_gene_enrich$gene)


Fig_3d <-  ggplot(SFARI_gene_enrich, aes(x = cluster, y = -log10(p_adjusted)) ) +
  geom_col(aes(x = cluster, y = -log10(p_adjusted), fill = cluster), 
                   colour = "black",
                     position = position_dodge()) +
  labs(y = expression(paste(-log[10], "(adj. ", italic("P"), ")")), title = Fig3d_title) +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept=-log10(0.05), linetype='dashed', col = "black", linewidth = 0.7) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values=safe_colorblind_palette) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        panel.grid.major.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 8, face = "bold"), 
        legend.position = "none" ) 
Fig_3d


#Export figure
Figure_3 <- ggdraw() +
  draw_plot(Fig_3A, x = 0, y = 0.63, width = 0.5, height = 0.375) +
  draw_plot(genesum_river_updated, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot(combined_plot, x = 0, y = 0.19, width = 0.5, height = 0.44) +
  draw_plot(Fig_3d, x = 0.0, y = 0, width = 0.5, height = 0.191) +
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 14, 
                  x = c(0, 0.5, 0, 0),
                  y = c(1, 1, 0.63, 0.1911))
#Figure_3

pdf(paste0("./figures/Figure_3.pdf"), height = 6.5, width = 6.5)
print(Figure_3)
dev.off()


svg(paste0("./figures/Figure_3.svg"), height = 6.5, width = 6.5)
print(Figure_3)
dev.off()


# FIGURE 4 --------------------------------------------------------------
classification <- read_parquet("./data/final_classification.parquet") #182371 

###IsoformSwitches
aSwitchList_part2 <- readRDS("./code/IsoformSwitchAnalyzeR/output/isoformswitch_part2_v4.rds")
aSwitchList_part2[["isoformFeatures"]]$gene_name <- aSwitchList_part2[["isoformFeatures"]]$gene_id




# FIGURE 4A ---------------------------------------------------------------
###DTU PROPORTIONS (bars)
DTU <- read.csv("./code/IsoformSwitchAnalyzeR/output/DTU_table.csv")
DTU_filter <- DTU %>% filter(DTU_qval < 0.05 & 
                               (DTU_dIF < -0.05 | DTU_dIF > 0.05)) 

DTU_filter <- DTU_filter[,-c(1:3)]
DTU_filter$condition <- paste0(DTU_filter$condition_1, " vs ", DTU_filter$condition_2)
DTU_filter <- merge(DTU_filter, classification[, c("isoform", "structural_category")],
     by.x = "isoform_id", by.y = "isoform", all.x = T)

rep_str = c('full-splice_match'='FSM','incomplete-splice_match' = 'ISM', 'novel_in_catalog'='NIC','novel_not_in_catalog'='NNC')
DTU_filter$structural_category <- str_replace_all(DTU_filter$structural_category, rep_str)

#All conditions
DTU_summ <- DTU_filter %>%
  group_by(condition, structural_category) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(perc = (count/sum(count))*100) %>%
  ungroup()

DTU_summ$percentage = paste0(round(DTU_summ$perc, 0), "%")

####DTE
res_t30_vs_t00 <- read.csv("./code/expression/output/DESeq2_tr_t30_vs_t00.csv")
res_t30_vs_t00$condition <- "t00 vs t30"

res_t04_vs_t00 <- read.csv("./code/expression/output/DESeq2_tr_t04_vs_t00.csv")
res_t04_vs_t00$condition <- "t00 vs t04"


res_t30_vs_t04 <- read.csv("./code/expression/output/DESeq2_tr_t30_vs_t04.csv")
res_t30_vs_t04$condition <- "t04 vs t30"


DTE_prefilter <- rbind(res_t30_vs_t00, res_t04_vs_t00)
DTE_prefilter <- rbind(DTE_prefilter, res_t30_vs_t04)


DTE_prefilter <- merge(DTE_prefilter, classification[, c("isoform", "associated_gene", "structural_category")],
                        by.x = "pb_id", by.y = "isoform", all.x = T)

#Just take the transcripts that meet the DTU pre-filter:
prefilter <- aSwitchList_part2$isoformRepIF$isoform_id

DTE_prefilter <- DTE_prefilter %>% filter(pb_id %in% prefilter) #68,801
DTE <- DTE_prefilter %>% filter(padj < 0.05 & 
                                   (log2FoldChange < -1 | log2FoldChange > 1))


DTE$structural_category <- str_replace_all(DTE$structural_category, rep_str)

DTE_summ <- DTE %>%
  group_by(condition, structural_category) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(perc = (count/sum(count))*100) %>%
  ungroup()

DTE_summ$percentage = paste0(round(DTE_summ$perc, 0), "%")


DTE_summ$type <- "DTE"
DTU_summ$type <- "DTU"
DTU_summ <- rbind(DTU_summ, DTE_summ)

              #FSM        NIC        NNC      ISM
colours <- c("#009E73", "#D55E00", "#E69F00", "#0072B2")



DTU_summ$condition <- factor(DTU_summ$condition, levels = c("t00 vs t04", "t04 vs t30", "t00 vs t30"))


bars_by_tr <- ggplot(DTU_summ, aes(x = condition, y = count/(1000), fill = structural_category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 1.9) +
  theme_classic(base_size = 8) +
  scale_fill_manual(values = colours[c(1,4,2,3)]) +
  theme(
    legend.title = element_blank(),
    legend.key.size = unit(0.35, 'cm'), 
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 6.4, colour = "black"), 
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"), 
    axis.text.x=element_text(angle=-45, hjust = 0, vjust=1)) + 
  labs(y = expression(paste("Isoforms (x ", 10^{3}, ")"))) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  facet_wrap(~factor(type), scales = "free")
bars_by_tr
#


# FUNC. CONSEQ ENRICHMENT------------------------------------------------------
#Need to run the function fully to extract genes, etc.
switchAnalyzeRlist = aSwitchList_part2
consequencesToAnalyze = 'all'
alpha=0.05
dIFcutoff = 0.05
countGenes = TRUE
analysisOppositeConsequence=FALSE
plot=TRUE
localTheme = theme_bw(base_size = 12)
minEventsForPlotting = 10
returnResult=TRUE
returnSummary=TRUE

### Consequences to analyze
acceptedTypes <- c(
      # Transcript
      'tss',
      'tts',
      'last_exon',
      'isoform_length',
      'exon_number',
      'intron_structure',
      'intron_retention',
      'isoform_class_code',
      # cpat
      'coding_potential',
      # ORF
      'ORF_genomic',
      'ORF_length',
      '5_utr_length',
      '3_utr_length',
      # seq similarity
      'isoform_seq_similarity',
      'ORF_seq_similarity',
      '5_utr_seq_similarity',
      '3_utr_seq_similarity',
      # ORF
      'NMD_status',
      # pfam
      'domains_identified',
      'genomic_domain_position',
      'domain_length',
      'domain_isotype',
      
      # SignalIP
      'signal_peptide_identified',
      
      # IDR
      'IDR_identified',
      'IDR_length',
      'IDR_type',
      
      # sub cell
      'sub_cell_location',
      'sub_cell_shift_to_cell_membrane',
      'sub_cell_shift_to_cytoplasm',
      'sub_cell_shift_to_nucleus',
      'sub_cell_shift_to_Extracellular'
      
      # topology
    #  'isoform_topology'
    )
    
consequencesAnalyzed <- unique(switchAnalyzeRlist$switchConsequence$featureCompared)
if ('all' %in% consequencesToAnalyze) {
      consequencesToAnalyze <- consequencesAnalyzed
    }
    
consequencesNotAnalyzed <- setdiff(consequencesToAnalyze, consequencesAnalyzed)
if (length(consequencesNotAnalyzed)) {
      warning(
        paste(
          'The following consequences appear not to have been analyzed and will therefor not be summarized:',
          paste(consequencesNotAnalyzed, collapse = ', '),
          sep = ' '
        )
      )
    }
  
### Extract non-location consequences
  if(TRUE) {
    localConseq <- switchAnalyzeRlist$switchConsequence[
      which( !is.na(
        switchAnalyzeRlist$switchConsequence$switchConsequence
      ))
      ,]
    localConseq <- localConseq[which(
      ! grepl('switch', localConseq$switchConsequence)
    ),]
    
    ### make list with levels
    levelList <- list(
      tss=c('Tss more upstream','Tss more downstream'),
      tts=c('Tts more downstream','Tts more upstream'),
      last_exon=c('Last exon more downstream','Last exon more upstream'),
      isoform_length=c('Length gain','Length loss'),
      isoform_seq_similarity=c('Length gain','Length loss'),
      exon_number=c('Exon gain','Exon loss'),
      intron_retention=c('Intron retention gain','Intron retention loss'),
      ORF_length=c('ORF is longer','ORF is shorter'),
      ORF=c('Complete ORF loss','Complete ORF gain'),
      x5_utr_length=c('5UTR is longer','5UTR is shorter'),
      x3_utr_length=c('3UTR is longer','3UTR is shorter'),
      NMD_status=c('NMD sensitive','NMD insensitive'),
      coding_potential=c('Transcript is coding','Transcript is Noncoding'),
      domains_identified=c('Domain gain','Domain loss'),
      domain_length=c('Domain length gain','Domain length loss'),
      domain_isotype=c('Domain non-reference isotype gain','Domain non-reference isotype loss'),
      IDR_identified = c('IDR gain','IDR loss'),
      IDR_length = c('IDR length gain', 'IDR length loss'),
      IDR_type = c('IDR w binding region gain', 'IDR w binding region loss'),
      signal_peptide_identified=c('Signal peptide gain','Signal peptide loss'),
      sub_cell_location = c('SubCell location gain','SubCell location loss'),
      sub_cell_shift_to_cell_membrane = c('SubCell location memb gain','SubCell location memb loss'),
      sub_cell_shift_to_cytoplasm     = c('SubCell location cyto gain','SubCell location cyto loss'),
      sub_cell_shift_to_nucleus       = c('SubCell location nucl gain','SubCell location nucl loss'),
      sub_cell_shift_to_Extracellular = c('SubCell location ext cell gain','SubCell location ext cell loss'),
      isoform_topology = c('Topology complexity gain','Topology complexity loss'),
      extracellular_region_count = c('Extracellular region gain', 'Extracellular region loss'),
      intracellular_region_count = c('Intracellular region gain', 'Intracellular region loss'),
      extracellular_region_length = c('Extracellular length gain','Extracellular length loss'),
      intracellular_region_length = c('Intracellular length gain','Intracellular length loss')
    )
    levelListDf <- plyr::ldply(levelList, function(x) data.frame(feature=x, stringsAsFactors = FALSE))
    
### Add consequence pairs
    localConseq$conseqPair <- levelListDf$.id[match(localConseq$switchConsequence, levelListDf$feature)]
    localConseq <- localConseq[which( !is.na(localConseq$conseqPair)),]
    
    ### Subset to consequences analyzed
    localConseq <- localConseq[which(
      localConseq$featureCompared %in% consequencesToAnalyze
    ),]
  }
  
### Subset to significant features
  if(TRUE) {
    ### Extract Sig iso
    isoResTest <-
      any(!is.na(
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
      ))
    if (isoResTest) {
      sigIso <- switchAnalyzeRlist$isoformFeatures[which(
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value < alpha &
          abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
      ),
      c('iso_ref', 'gene_ref')]
    } else {
      sigIso <- switchAnalyzeRlist$isoformFeatures[which(
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
          abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
      ),
      c('iso_ref', 'gene_ref')]
    }
    
    if(isoResTest) {
      localConseq <- localConseq[which(
        localConseq$iso_ref_down %in% sigIso$iso_ref |
          localConseq$iso_ref_up   %in% sigIso$iso_ref
      ),]
    } else {
      localConseq <- localConseq[which(
        localConseq$gene_ref %in% sigIso$gene_ref
      ),]
    }
  }
  
  
### Summarize gain vs loss for each consequence in each condition
  consequenceBalance <- plyr::ddply(
    .data = localConseq,
    .variables = c('condition_1','condition_2','conseqPair'),
    #.inform = TRUE,
    .fun = function(aDF) { # aDF <- localConseq[1:20,]
      ### Add levels
      if(analysisOppositeConsequence) {
        localLvl <- rev(sort(
          levelList[[ aDF$conseqPair[1] ]]
        ))
      } else {
        localLvl <- sort(
          levelList[[ aDF$conseqPair[1] ]]
        )
      }
      aDF$switchConsequence <- factor(
        aDF$switchConsequence,
        levels=localLvl
      )
      
      ### Summarize category
      if( countGenes ) {
        df2 <- aDF[
          which(!is.na(aDF$switchConsequence)),
          c('gene_id','switchConsequence')
        ]
        localNumber <- plyr::ddply(df2, .drop = FALSE, .variables = 'switchConsequence', function(x) {
          data.frame(
            Freq = length(unique(x$gene_id))
          )
        })
        colnames(localNumber)[1] <- 'Var1'
      } else {
        localNumber <- as.data.frame(table(aDF$switchConsequence))
      }
      
      if(nrow(localNumber) == 2) {
        localTest <- suppressWarnings(
          stats::binom.test(localNumber$Freq[1], sum(localNumber$Freq))
        )
        
        localRes <- data.frame(
          feature=stringr::str_c(
            localNumber$Var1[1],
            ' (paired with ',
            localNumber$Var1[2],
            ')'
          ),
          propOfRelevantEvents=localTest$estimate,
          stringsAsFactors = FALSE
        )
        
        localRes$propCiLo <- min(localTest$conf.int)
        localRes$propCiHi <- max(localTest$conf.int)
        localRes$propPval <- localTest$p.value
      } else {
        warning('Somthing strange happend - contact developer with reproducible example')
      }
      
      localRes$nUp   <- localNumber$Freq[which( localNumber$Var1 == levels(localNumber$Var1)[1] )]
      localRes$nDown <- localNumber$Freq[which( localNumber$Var1 == levels(localNumber$Var1)[2] )]
      
      return(localRes)
    }
  )
  
  consequenceBalance$propQval <- p.adjust(consequenceBalance$propPval, method = 'fdr')
  consequenceBalance$Significant <- consequenceBalance$propQval < alpha
  consequenceBalance$Significant <- factor(
    consequenceBalance$Significant,
    levels=c(FALSE,TRUE)
  )

#Prep for plotting
  consequenceBalance2 <- consequenceBalance[which(
      (consequenceBalance$nUp + consequenceBalance$nDown) >= minEventsForPlotting
    ),]
    if(nrow(consequenceBalance2) == 0) {
      stop('No features left for ploting after filtering with via "minEventsForPlotting" argument.')
    }

    consequenceBalance2$nTot <- consequenceBalance2$nDown + consequenceBalance2$nUp

    ### Add comparison
    consequenceBalance2$Comparison <- paste(
      consequenceBalance2$condition_1,
      'vs',
      consequenceBalance2$condition_2,
      sep='\n'
    )

    ### Massage
    consequenceBalance2$feature2 <- gsub(' \\(', '\n(', consequenceBalance2$feature)

    consequenceBalance2$feature2 <- factor(
      consequenceBalance2$feature2,
      levels = rev(sort(unique(as.character(consequenceBalance2$feature2))))
    )

# SPLICING PLOT -----------------------------------------------------------
extractSwitchPairs <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.05,
    onlySigIsoforms = FALSE
  ) {
    ### Extract and massage data
    if (TRUE) {
      localData <- switchAnalyzeRlist$isoformFeatures[
        which(
          switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
            abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
        ),
        c(
          'iso_ref',
          'gene_ref',
          'isoform_switch_q_value',
          'gene_switch_q_value',
          'dIF'
        )
      ]
      
      if (!nrow(localData)) {
        stop('No genes were considered switching with the used cutoff values')
      }
      
      ### add switch direction
      localData$switchDirection <- NA
      localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
      localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'
      
      ### Annotate significant features
      isoResTest <-
        any(!is.na(
          switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
        ))
      if (isoResTest) {
        localData$isoSig <-
          localData$isoform_switch_q_value < alpha &
          abs(localData$dIF) > dIFcutoff
      } else {
        localData$isoSig <-
          localData$gene_switch_q_value < alpha &
          abs(localData$dIF) > dIFcutoff
      }
      
      
      if(onlySigIsoforms) {
        localData <- localData[which( localData$isoSig ),]
      }
    }
    
    
    ### Create data sub-sets of interest
    if(TRUE) {
      sigUpData <- localData[which(
        localData$isoSig & localData$switchDirection == 'up'
      ),c('iso_ref','gene_ref')]
      sigDnData <- localData[which(
        localData$isoSig & localData$switchDirection == 'down'
      ),c('iso_ref','gene_ref')]
      
      colnames(sigUpData)[1] <- c('iso_ref_up')
      colnames(sigDnData)[1] <- c('iso_ref_down')
      
      
      if( ! onlySigIsoforms ) {
        justUpData <- localData[which(
          localData$switchDirection == 'up'
        ),c('iso_ref','gene_ref')]
        justDnData <- localData[which(
          localData$switchDirection == 'down'
        ),c('iso_ref','gene_ref')]
        
        colnames(justUpData)[1]  <- c('iso_ref_up')
        colnames(justDnData)[1] <- c('iso_ref_down')
      }
    }
    
    ### Join datasets to extract pairs
    if(TRUE) {
      if( onlySigIsoforms ) {
        pairwiseIsoComparison <- dplyr::inner_join(
          sigUpData,
          sigDnData,
          by= 'gene_ref',
          multiple = "all"
        )
      } else {
        ### Sig up and all down
        upPairs <- dplyr::inner_join(
          sigUpData,
          justDnData,
          by= 'gene_ref',
          multiple = "all"
        )
        ### Sig down and all up
        dnPairs <- dplyr::inner_join(
          justUpData,
          sigDnData,
          by= 'gene_ref',
          multiple = "all"
        )
        
        ### Combine
        pairwiseIsoComparison <- unique(
          rbind(
            upPairs,
            dnPairs
          )
        )
        pairwiseIsoComparison <- pairwiseIsoComparison[,c(
          'gene_ref','iso_ref_up','iso_ref_down'
        )]
        
        ### Reorder
        pairwiseIsoComparison <- pairwiseIsoComparison[order(
          pairwiseIsoComparison$gene_ref,
          pairwiseIsoComparison$iso_ref_up,
          pairwiseIsoComparison$iso_ref_down
        ),]
      }
    }
    
    ### Add in additional data
    if(TRUE) {
      ### Add isoform names
      matchVectorUp <- match(
        pairwiseIsoComparison$iso_ref_up,
        switchAnalyzeRlist$isoformFeatures$iso_ref
      )
      matchVectorDn <- match(
        pairwiseIsoComparison$iso_ref_down,
        switchAnalyzeRlist$isoformFeatures$iso_ref
      )
      
      pairwiseIsoComparison$isoformUpregulated   <-
        switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorUp]
      pairwiseIsoComparison$isoformDownregulated <-
        switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorDn]
      
      ### Gene infl
      pairwiseIsoComparison$gene_id <-
        switchAnalyzeRlist$isoformFeatures$gene_id[matchVectorUp]
      pairwiseIsoComparison$gene_name <-
        switchAnalyzeRlist$isoformFeatures$gene_name[matchVectorUp]
      
      ### Conditons
      pairwiseIsoComparison$condition_1 <-
        switchAnalyzeRlist$isoformFeatures$condition_1[matchVectorUp]
      pairwiseIsoComparison$condition_2 <-
        switchAnalyzeRlist$isoformFeatures$condition_2[matchVectorUp]
      
    }
    
    return( pairwiseIsoComparison )
  }
  
#Use MIC IsoformSwitch object to extract the MICs
aSwitchList_part2_MIC <- readRDS(file="./code/microexons/output/isoformswitch_MIC.rds") 

switchAnalyzeRlist = aSwitchList_part2_MIC
splicingToAnalyze = 'all'
alpha = 0.05
dIFcutoff = 0.05
onlySigIsoforms = T
countGenes = T
plot = TRUE
localTheme = theme_bw(base_size = 14)
minEventsForPlotting = 25 ##default = 10
returnResult=TRUE
returnSummary=TRUE

### Consequences to analyze
acceptedTypes <- c("A3","A5","ATSS","ATTS","ES" ,"IR","MEE","MES", "MIC" )


splicingAnalyzed <-
  intersect(
    acceptedTypes,
    colnames(switchAnalyzeRlist$AlternativeSplicingAnalysis)
  )
if ('all' %in% splicingToAnalyze) {
  splicingToAnalyze <- splicingAnalyzed
}

splicingNotAnalyzed <-
  setdiff(splicingToAnalyze, splicingAnalyzed)
if (length(splicingNotAnalyzed)) {
  warning(
    paste(
      'The following consequences appear not to have been analyzed and will therefor not be summarized:',
      paste(splicingNotAnalyzed, collapse = ', '),
      sep = ' '
    )
  )
}

### Get pairs
if(TRUE) {
  pairwiseIsoComparison <- extractSwitchPairs(
    switchAnalyzeRlist,
    alpha = alpha,
    dIFcutoff = dIFcutoff,
    onlySigIsoforms = onlySigIsoforms
  )
}

### Massage AS analysis
if(TRUE) {
  localAS <- switchAnalyzeRlist$AlternativeSplicingAnalysis
  localAS <- localAS[which(
    localAS$isoform_id %in% pairwiseIsoComparison$isoformUpregulated |
      localAS$isoform_id %in% pairwiseIsoComparison$isoformDownregulated
  ),]
  
  ### Massage
  m1 <- reshape2::melt(localAS[,c(
    "isoform_id",
    "ES_genomic_start",
    "MEE_genomic_start",
    "MES_genomic_start",
    "IR_genomic_start",
    "A5_genomic_start",
    "A3_genomic_start",
    "ATSS_genomic_start",
    "ATTS_genomic_start", 
    "MIC_genomic_start"
  )], id.vars = 'isoform_id')
  colnames(m1)[3] <- 'genomic_start'
  m1$AStype <- sapply(
    strsplit(as.character(m1$variable),'_'),
    function(x) x[1]
  )
  
  m2 <- reshape2::melt(localAS[,c(
    "isoform_id",
    "ES_genomic_end",
    "MEE_genomic_end",
    "MES_genomic_end",
    "IR_genomic_end",
    "A5_genomic_end",
    "A3_genomic_end",
    "ATSS_genomic_end",
    "ATTS_genomic_end",
    "MIC_genomic_end"
  )], id.vars = 'isoform_id')
  colnames(m2)[3] <- 'genomic_end'
  m2$AStype <- sapply(
    strsplit(as.character(m2$variable),'_'),
    function(x) x[1]
  )
  
  localAS <- dplyr::inner_join(
    m1[,c('isoform_id','AStype','genomic_start')],
    m2[,c('isoform_id','AStype','genomic_end')],
    by=c('isoform_id','AStype')
  )
  
  ### Add in NMD
  if('orfAnalysis' %in% names(switchAnalyzeRlist)) {
    localNMD <- data.frame(
      isoform_id=switchAnalyzeRlist$orfAnalysis$isoform_id,
      AStype='NMD',
      genomic_start=ifelse(switchAnalyzeRlist$orfAnalysis$PTC, '0,0', '0'),
      genomic_end=ifelse(switchAnalyzeRlist$orfAnalysis$PTC, '0,0', '0'),
      stringsAsFactors = FALSE
    )
    localAS <- rbind(localAS , localNMD)
  }
}

### Add AS to pairs
if(TRUE) {
  localConseq2 <- merge(
    pairwiseIsoComparison,
    localAS,
    by.x='isoformUpregulated',
    by.y='isoform_id'
  )
  
  localConseq3 <- merge(
    localConseq2,
    localAS,
    by.x=c('isoformDownregulated','AStype'),
    by.y=c('isoform_id','AStype'),
    suffixes = c("_up","_down")
  )
  
  
  ### Replace NAs so they can be compared (na just mean same as pre-transcript)
  localConseq3$genomic_start_up[which(
    is.na(localConseq3$genomic_start_up)
  )] <- 0
  localConseq3$genomic_end_up[which(
    is.na(localConseq3$genomic_end_up)
  )] <- 0
  localConseq3$genomic_start_down[which(
    is.na(localConseq3$genomic_start_down)
  )] <- 0
  localConseq3$genomic_end_down[which(
    is.na(localConseq3$genomic_end_down)
  )] <- 0
  
  ### Identify differences
  localConseq3$coordinatsDifferent <-
    localConseq3$genomic_start_up != localConseq3$genomic_start_down |
    localConseq3$genomic_end_up   != localConseq3$genomic_end_down
  
  localConseq4 <- localConseq3[which(localConseq3$coordinatsDifferent),]
  
  if(nrow(localConseq4) == 0) {
    stop('No alternative splicing differences were found')
  }
  
}

### Subset based on input
if(TRUE) {
  localConseq4 <- localConseq4[which(
    localConseq4$AStype %in% splicingToAnalyze
  ),]
  if (!nrow(localConseq4)) {
    stop('No swithces with consequences were found')
  }
  
}

### Analyze differences (both start and end coordinats)
if(TRUE) {
  genomic_start_up   <- strsplit(x = localConseq4$genomic_start_up  , split = ';')
  genomic_start_down <- strsplit(x = localConseq4$genomic_start_down, split = ';')
  genomic_end_up   <- strsplit(x = localConseq4$genomic_end_up  , split = ';')
  genomic_end_down <- strsplit(x = localConseq4$genomic_end_down, split = ';')
  
  localConseq4$nrGain <- 0
  localConseq4$nrLoss <- 0
  for(i in seq_along(genomic_start_up)) { 
    # combine start and end of exons
    localUp <- stringr::str_c(genomic_start_up[[i]]  , '_', genomic_end_up  [[i]])
    localDn <- stringr::str_c(genomic_start_down[[i]], '_', genomic_end_down[[i]])
    
    # remove coordinats from primary transcripts
    localUp <- localUp[which( localUp != '0_0')]
    localDn <- localDn[which( localDn != '0_0')]
    
    # calculate difference
    localConseq4$nrGain[i] <- sum(!localUp %in% localDn)
    localConseq4$nrLoss[i] <- sum(!localDn %in% localUp)
  }
  ### Modify ends (can only have 1)
  terminiIndex <- which( localConseq4$AStype %in% c('ATSS','ATTS') )
  localConseq4$nrGain[terminiIndex] <-
    1 * sign(localConseq4$nrGain[terminiIndex])
  localConseq4$nrLoss[terminiIndex] <-
    1 * sign(localConseq4$nrLoss[terminiIndex])
  
  localConseq4$nrDiff <- localConseq4$nrGain - localConseq4$nrLoss
  
}

### Summarize gain vs loss for each AStype in each condition
gainLossBalance <- plyr::ddply(
  .data = localConseq4,
  .variables = c('condition_1','condition_2','AStype'),
  .fun = function(aDF) { # aDF <- localConseq4[1:50,]
    if( countGenes ) {
      aDF2 <- plyr::ddply(aDF[,c('gene_ref','nrDiff')], .variables = 'gene_ref', function(aGene) {
        data.frame(
          nrDiff = sum(aGene$nrDiff)
        )
      })
      
      localRes <- data.frame(
        nUp  = sum(aDF2$nrDiff > 0),
        nDown= sum(aDF2$nrDiff < 0)
      )
    } else {
      localRes <- data.frame(
        nUp=sum(aDF$nrDiff > 0),
        nDown=sum(aDF$nrDiff < 0)
      )
    }
    
    if( localRes$nUp > 0 | localRes$nDown > 0 ) {
      localTest <- suppressWarnings(
        stats::binom.test(localRes$nUp, localRes$nUp + localRes$nDown)
      )
      
      localRes$propUp <- localTest$estimate
      localRes$propUpCiLo <- min(localTest$conf.int)
      localRes$propUpCiHi <- max(localTest$conf.int)
      localRes$propUpPval <- localTest$p.value
    } else {
      localRes$propUp <- NA
      localRes$propUpCiLo <- NA
      localRes$propUpCiHi <- NA
      localRes$propUpPval <- NA
    }
    
    return(localRes)
  }
)

### Massage for plotting
if(TRUE) {
  gainLossBalance <- gainLossBalance[which(
    ! is.na(gainLossBalance$propUp)
  ),]
  gainLossBalance$propUpQval <- p.adjust(gainLossBalance$propUpPval, method = 'fdr')
  gainLossBalance$Significant <- gainLossBalance$propUpQval < alpha
  gainLossBalance$Significant <- factor(
    gainLossBalance$Significant,
    levels=c(FALSE,TRUE)
  )
  
  gainLossBalance$Comparison <- paste(
    gainLossBalance$condition_1,
    'vs',
    gainLossBalance$condition_2,
    sep='\n'
  )
  
  ### Massage splicing type name
  if(TRUE) {
    gainLossBalance$AStype <- paste0(
      gainLossBalance$AStype, ' gain ',
      '(paired with ',gainLossBalance$AStype, ' loss)'
    )
    
    gainLossBalance$AStype <- gsub('MES gain', 'MES', gainLossBalance$AStype)
    gainLossBalance$AStype <- gsub('MES loss', 'MEI', gainLossBalance$AStype)
    
    gainLossBalance$AStype <- gsub('ES gain', 'ES', gainLossBalance$AStype)
    gainLossBalance$AStype <- gsub('ES loss', 'EI', gainLossBalance$AStype)
    
    #gainLossBalance$AStype <- gsub('IR gain', 'IR', gainLossBalance$AStype)
    #gainLossBalance$AStype <- gsub('IR loss', 'IS', gainLossBalance$AStype)
  }
  
  myOrder <- plyr::ddply(
    gainLossBalance,
    .variables = 'AStype',
    .fun = function(aDF) {
      mean(aDF$propUp)
    })
  myOrder <- myOrder$AStype[sort.list(myOrder$V1, decreasing = TRUE)]
  
  gainLossBalance$AStype <- factor(
    gainLossBalance$AStype,
    levels = myOrder
  )
}
# 
# Prep data for plotting
gainLossBalance2 <- gainLossBalance[which(
    (gainLossBalance$nUp + gainLossBalance$nDown) >= minEventsForPlotting
  ),]

  gainLossBalance2$nTot <- gainLossBalance2$nUp + gainLossBalance2$nDown


####Combining the data to plot:
opposite <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = TRUE,
  analysisOppositeConsequence=T,
  plot=F,
  returnResult=TRUE,
  returnSummary=TRUE)
opposite <- opposite %>% filter(conseqPair == "NMD_status" |
                                  conseqPair == "domains_identified" |
                                  conseqPair == "IDR_identified")
opposite$Conseq_type <- "Functional Consequence"
opposite$nTot <- opposite$nUp + opposite$nDown
opposite$Comparison <- paste0(opposite$condition_1, " vs ", opposite$condition_2)
opposite <- opposite %>% relocate(c(nUp, nDown), .after = feature)
opposite <- opposite %>% relocate(Comparison, .before = nTot)
opposite <- opposite[,-c(3)]
opposite <- opposite %>% relocate(c(Conseq_type), .after = nTot)


combined_data <- consequenceBalance2[,-c(3,15)]
combined_data$Conseq_type <- "Functional Consequence"
combined_data <- combined_data %>% relocate(c(nUp, nDown), .after = feature)
combined_data <- combined_data %>% relocate(Comparison, .before = nTot)
combined_data <- combined_data %>% filter(feature != "NMD insensitive (paired with NMD sensitive)" &
                                            feature != "Domain gain (paired with Domain loss)" &
                                            feature != "IDR gain (paired with IDR loss)")


colnames(combined_data) <- c("condition_1", "condition_2", "feature", "nUp", "nDown", 
                             "propUp" , "propUpCiLo", "propUpCiHi", "propUpPval", 
                             "propUpQval",  "Significant", "Comparison" , "nTot", "Conseq_type")
gainLossBalance2$Conseq_type <- "Alternative Splicing"
colnames(gainLossBalance2) <- c("condition_1", "condition_2", "feature", "nUp", "nDown", 
                             "propUp" , "propUpCiLo", "propUpCiHi", "propUpPval", 
                             "propUpQval",  "Significant", "Comparison" , "nTot", "Conseq_type")
colnames(opposite) <- c("condition_1", "condition_2", "feature", "nUp", "nDown", 
                             "propUp" , "propUpCiLo", "propUpCiHi", "propUpPval", 
                             "propUpQval",  "Significant", "Comparison" , "nTot", "Conseq_type")
combined_data <- rbind(combined_data, gainLossBalance2)
combined_data <- rbind(combined_data, opposite)

combined_data$Comparison <- str_replace(combined_data$Comparison, "t00\nvs\nt04", "t00 vs t04")
combined_data$Comparison <- str_replace(combined_data$Comparison, "t04\nvs\nt30", "t04 vs t30")
combined_data$Comparison <- str_replace(combined_data$Comparison, "t00\nvs\nt30", "t00 vs t30")

combined_data$feature <- as.factor(combined_data$feature)
library(forcats)


combined_data$feature <- str_replace(combined_data$feature, "3UTR is longer \\(paired with 3UTR is shorter\\)",
                                                            "3' UTR is longer")
combined_data$feature <- str_replace(combined_data$feature, "5UTR is longer \\(paired with 5UTR is shorter\\)", 
                                                            "5' UTR is longer")
combined_data$feature <- str_replace(combined_data$feature, "Domain loss \\(paired with Domain gain\\)",
                                                            "Domain loss")
combined_data$feature <- str_replace(combined_data$feature, "Domain length gain \\(paired with Domain length loss\\)",
                                                            "Domain length gain")
combined_data$feature <- str_replace(combined_data$feature, "Domain non-reference isotype gain \\(paired with Domain non-reference isotype loss\\)",
                                                            "Domain non-reference isotype gain")
combined_data$feature <- str_replace(combined_data$feature, "Exon gain \\(paired with Exon loss\\)", 
                                                            "Exon number gain")
combined_data$feature <- str_replace(combined_data$feature, "IDR loss \\(paired with IDR gain\\)",
                                                            "IDR loss")
combined_data$feature <- str_replace(combined_data$feature, "IDR length gain \\(paired with IDR length loss\\)",
                                                            "IDR is longer")
combined_data$feature <- str_replace(combined_data$feature, "IDR w binding region gain \\(paired with IDR w binding region loss\\)",
                                                            "IDR with binding region gain")
combined_data$feature <- str_replace(combined_data$feature, "Last exon more downstream \\(paired with Last exon more upstream\\)",
                                                            "Last exon more downstream")
combined_data$feature <- str_replace(combined_data$feature, "Length gain \\(paired with Length loss\\)",
                                                            "Isoform is longer")
combined_data$feature <- str_replace(combined_data$feature, "NMD sensitive \\(paired with NMD insensitive\\)",
                                                            "NMD sensitive")
combined_data$feature <- str_replace(combined_data$feature, "ORF is longer \\(paired with ORF is shorter\\)",
                                                            "ORF is longer")
combined_data$feature <- str_replace(combined_data$feature, "Tss more downstream \\(paired with Tss more upstream\\)",
                                                            "TSS more downstream")
combined_data$feature <- str_replace(combined_data$feature, "Tts more downstream \\(paired with Tts more upstream\\)",
                                                            "TTS more downstream")
combined_data$feature <- str_replace(combined_data$feature, "A3 gain \\(paired with A3 loss\\)",
                                                            "3' splice site more downstream")
combined_data$feature <- str_replace(combined_data$feature, "A5 gain \\(paired with A5 loss\\)",
                                                            "5' splice site more upstream")
combined_data$feature <- str_replace(combined_data$feature, "ATSS gain \\(paired with ATSS loss\\)",
                                                            "ATSS gain")
combined_data$feature <- str_replace(combined_data$feature, "ATTS gain \\(paired with ATTS loss\\)",
                                                            "ATTS gain")
combined_data$feature <- str_replace(combined_data$feature, "ES \\(paired with EI\\)",
                                                            "Exon skipping")
combined_data$feature <- str_replace(combined_data$feature, "MIC gain \\(paired with MIC loss\\)",
                                                            "Microexon skipping")
combined_data$feature <- str_replace(combined_data$feature, "IR gain \\(paired with IR loss\\)",
                                                            "Intron retention")
combined_data$feature <- str_replace(combined_data$feature, "MES \\(paired with MEI\\)",
                                                            "Multiple exon skipping")

combined_data$feature <- str_replace(combined_data$feature, "Signal peptide gain \\(paired with Signal peptide loss\\)",
                                     "Signal peptide gain")
combined_data$feature <- str_replace(combined_data$feature, "SubCell location cyto gain \\(paired with SubCell location cyto loss\\)",
                                     "Shift to cytoplasm localization")
combined_data$feature <- str_replace(combined_data$feature, "SubCell location ext cell gain \\(paired with SubCell location ext cell loss\\)",
                                     "Shift to extracellular localization")
combined_data$feature <- str_replace(combined_data$feature, "SubCell location gain \\(paired with SubCell location loss\\)",
                                     "Subcellular localization gain")
combined_data$feature <- str_replace(combined_data$feature, "SubCell location memb gain \\(paired with SubCell location memb loss\\)",
                                     "Shift to cell membrane localization")
combined_data$feature <- str_replace(combined_data$feature, "SubCell location nucl gain \\(paired with SubCell location nucl loss\\)",
                                     "Shift to nuclear localization")

combined_data$feature <- str_replace(combined_data$feature, "Topology complexity gain \\(paired with Topology complexity loss\\)",
                                     "Topology complexity gain")
combined_data$feature <- str_replace(combined_data$feature, "Extracellular length gain \\(paired with Extracellular length loss\\)",
                                     "Extracellular region length gain")
combined_data$feature <- str_replace(combined_data$feature, "Extracellular region gain \\(paired with Extracellular region loss\\)",
                                     "Extracellular region gain")
combined_data$feature <- str_replace(combined_data$feature, "Intracellular length gain \\(paired with Intracellular length loss\\)",
                                     "Intracellular region length gain")
combined_data$feature <- str_replace(combined_data$feature, "Intracellular region gain \\(paired with Intracellular region loss\\)",
                                     "Intracellular region gain")



combined_data <- combined_data %>%
  mutate(feature = fct_relevel(feature, 
                               "3' splice site more downstream", "5' splice site more upstream",
                               "ATSS gain", "ATTS gain", "Exon skipping", "Microexon skipping", "Intron retention", "Multiple exon skipping",
                              # "Intracellular region gain",
                              # "Extracellular region gain",
                              # "Intracellular region length gain",
                              # "Extracellular region length gain",
                              # "Topology complexity gain",
                               "Shift to extracellular localization",
                               "Shift to cell membrane localization",
                               "Shift to nuclear localization",
                               "Shift to cytoplasm localization",
                               "Subcellular localization gain",
                               "Signal peptide gain",
                               "IDR with binding region gain", 
                               "IDR is longer", 
                               "IDR loss", 
                               "Domain non-reference isotype gain", 
                               "Domain length gain", 
                               "Domain loss", 
                               "NMD sensitive", 
                               "ORF is longer", 
                               "Exon number gain", 
                               "Isoform is longer",
                               "Last exon more downstream", 
                               "3' UTR is longer", 
                               "TTS more downstream",
                               "5' UTR is longer", 
                               "TSS more downstream"))

colnames(combined_data)[13] <- 'Genes'

library(ggh4x)
library(extrafont)
library(scales)
library(grid)

font_import(paths = "C:/Windows/Fonts", pattern = "cambria", prompt = F)  # This imports all the fonts from your system
loadfonts(device = "win")  # For Windows
windowsFonts()

plot_data <- combined_data %>% filter(#feature == "Multiple exon skipping" |
                           feature == "Intron retention" |
                           feature == "Exon skipping" |
                           feature == "Microexon skipping" |
                         #  feature == "5' splice site more upstream" |
                         #  feature == "3' splice site more downstream" |
                         #  feature == "TSS more downstream" |
                           feature == "5' UTR is longer" |
                         #  feature == "TTS more downstream" |
                           feature == "3' UTR is longer" |
                      #     feature == "Isoform is longer" |
                         #  feature == "ORF is longer" |
                           feature == "NMD sensitive" |
                           feature == "Domain loss" |
                      #     feature == "Domain length gain" |
                           feature == "IDR loss"  #|
                      #     feature == "IDR is longer"
                      )
desired_order <- c("Exon skipping",
                   "Microexon skipping",
                   "Intron retention",
                   "NMD sensitive",
                   "Domain loss", 
                   "IDR loss", 
                   "5' UTR is longer",
                   "3' UTR is longer")
plot_data <- plot_data %>%
  mutate(feature = factor(feature, levels = rev(desired_order)))


combined_plot_1 <- plot_data %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  # geom_point(aes(size=Genes)) +
  geom_text(label = "\u25D6", aes(size=Genes), family = "Cambria") +
  facet_wrap(~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_x") +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x="Fraction of Genes with Switches Primarily Resulting in Each Event") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.text = element_text(size = 7),
        legend.key.spacing.y = unit(0, "pt"),
        legend.key.height = unit(-20, "lines"),
        axis.title.y = element_blank(), 
        legend.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "cm"
        )) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1, max(plot_data$Genes)), 
               breaks = c(50, 500, 2000, 4000)) +
  guides(color = guide_legend(order=1, nrow = 2, byrow = T),
         size = guide_legend(order=2, override.aes = list(color = "black"), nrow = 2, byrow = T)) +
  coord_cartesian(xlim=c(0,1)) 

combined_plot_1

# Create a version of the plot for legend extraction
combined_plot_symbols_only <- combined_plot_1 +
  theme(
    legend.text = element_text(colour = "transparent"),  # invisible but still occupies space
    legend.title = element_text(colour = "transparent")  # invisible but still occupies space
  )

legend_1 <- ggpubr::get_legend(combined_plot_symbols_only)



combined_plot_2 <- plot_data %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  geom_text(label = "\u25D7", aes(size=Genes), family = "Cambria") +
  facet_wrap(~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_x") +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x="Fraction of Genes with Switches Primarily Resulting in Each Event"
  ) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.key.height = unit(-10, "lines"),
        legend.key.spacing.y = unit(0, "pt"),
        legend.text = element_text(size = 7),
        axis.title.y = element_blank(), 
        legend.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "cm"
        )) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1, max((plot_data$Genes))), 
               breaks = c(50, 500, 2000, 4000)) +
  guides(color = guide_legend(order=1, nrow = 2, byrow = T),
         size = guide_legend(order=2, override.aes = list(color = "red"), nrow = 2, byrow = T)) +
  coord_cartesian(xlim=c(0,1)) 

combined_plot_2


# Extract the legend from the plot
legend_2 <- ggpubr::get_legend(combined_plot_2)

combined_plot <- plot_data %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  # geom_point(aes(size=Genes)) +
  geom_text(label = "\u25CF", aes(size=Genes), family = "Cambria") +
  facet_wrap(~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_x") +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x=str_wrap("Fraction of Genes with Switches Primarily Resulting in Each Event", 40)
       # y='Event/Consequence of Isoform Switch'
  ) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 17, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "none",
        axis.title.y = element_blank()  ) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1, max((plot_data$Genes)))) +
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2)) +
  coord_cartesian(xlim=c(0,1)) +
  scale_y_discrete(labels = label_wrap_gen(width = 17)) 

combined_plot

testp <- ggdraw() +
  draw_plot(combined_plot + theme(plot.margin = margin(t = 15, r = 6, b = 6, l = 10, unit = "pt"))
            , x = 0, y = 0, width = 1, height = 0.94) +
  draw_plot(legend_1, x = 0, y = 0.41, width = 1, height = 1) +
  draw_plot(legend_2, x = 0, y = 0.41, width = 1, height = 1)
testp




# PLOTTING SWITCH SUPP FIG ------------------------------------------------
plot_data_full <- combined_data %>% filter(feature != "Shift to extracellular localization" &
    feature != "Shift to cell membrane localization" &
    feature != "Shift to nuclear localization" &
    feature != "Shift to cytoplasm localization" &
    feature != "Subcellular localization gain" &
    feature != "Signal peptide gain" &
    feature != "MEE gain (paired with MEE loss)"  )


combined_plot_1 <- plot_data_full %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  # geom_point(aes(size=Genes)) +
  geom_text(label = "\u25D6", aes(size=Genes), family = "Cambria") +
  facet_grid(Conseq_type ~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_y") +
  force_panelsizes(c((plot_data %>% filter(Conseq_type == "Alternative Splicing") %>% tally() /3),
                     (plot_data %>% filter(Conseq_type == "Functional Consequence") %>% tally() /3)), 
                   cols = 3, 
                   respect = NULL) +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x="Fraction of Genes with Switches Primarily Resulting in Each Event") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.text = element_text(size = 7),
        #  legend.key.size = element_rect(size = 10), 
        legend.key.spacing.y = unit(0, "pt"),
        legend.key.height = unit(-20, "lines"),
        axis.title.y = element_blank(), 
        legend.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "cm"
        )) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1,    6501), 
               breaks = c(50, 500, 2500, 6500)) +
  guides(color = guide_legend(order=1, nrow = 2, byrow = T),
         size = guide_legend(order=2, override.aes = list(color = "black"), nrow = 2, byrow = T)) +
  coord_cartesian(xlim=c(0,1)) 

combined_plot_1


# Create a version of the plot for legend extraction
combined_plot_symbols_only <- combined_plot_1 +
  theme(
    legend.text = element_text(colour = "transparent"),  # invisible but still occupies space
    legend.title = element_text(colour = "transparent")  # invisible but still occupies space
  )

legend_1 <- ggpubr::get_legend(combined_plot_symbols_only)



combined_plot_2 <- plot_data_full %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  geom_text(label = "\u25D7", aes(size=Genes), family = "Cambria") +
  facet_grid(Conseq_type ~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_y") +
  force_panelsizes(c((plot_data %>% filter(Conseq_type == "Alternative Splicing") %>% tally() /3),
                     (plot_data %>% filter(Conseq_type == "Functional Consequence") %>% tally() /3)), #7 AS, 15 func conseq
                   cols = 3, 
                   respect = NULL) +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x="Fraction of Genes with Switches Primarily Resulting in Each Event"
       # y='Event/Consequence of Isoform Switch'
  ) +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "top",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.key.height = unit(-10, "lines"),
        legend.key.spacing.y = unit(0, "pt"),
        legend.text = element_text(size = 7),
        axis.title.y = element_blank(), 
        legend.margin = margin(
          t = 0,
          r = 0,
          b = 0,
          l = 0,
          unit = "cm"
        )) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1,   6501), 
               breaks = c(50, 500, 2500, 6500)) +
  guides(color = guide_legend(order=1, nrow = 2, byrow = T),
         size = guide_legend(order=2, override.aes = list(color = "red"), nrow = 2, byrow = T)) +
  coord_cartesian(xlim=c(0,1)) 

combined_plot_2


# Extract the legend from the plot
legend_2 <- ggpubr::get_legend(combined_plot_2)

combined_plot <- plot_data_full %>%
  ggplot(aes(y=feature, x=propUp, color=Significant)) +
  geom_errorbarh(aes(xmax = propUpCiLo, xmin=propUpCiHi), height = .3) +
  # geom_point(aes(size=Genes)) +
  geom_text(label = "\u25CF", aes(size=Genes), family = "Cambria") +
  facet_grid(Conseq_type ~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")), scales = "free_y") +
  force_panelsizes(c((plot_data %>% filter(Conseq_type == "Alternative Splicing") %>% tally() /3),
                     (plot_data %>% filter(Conseq_type == "Functional Consequence") %>% tally() /3)), #7 AS, 15 func conseq
                   cols = 3, 
                   respect = NULL) +
  geom_vline(xintercept=0.5, linetype='dashed') +
  labs(title = "Isoform Switching Consequences",
       x=str_wrap("Fraction of Genes with Switches Primarily Resulting in Each Event", 40)
       # y='Event/Consequence of Isoform Switch'
  ) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 17, size = 8, face = "bold"),
        strip.text = element_text(size = 8),
        axis.text = element_text(colour = "black", size = 6.4),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust=1), 
        legend.position = "none",
        axis.title.y = element_blank()  ) +
  scale_color_manual(name = paste0('FDR < ', alpha), values=c('black','red'), drop=FALSE) +
  scale_radius(limits=c(1, 6501)) +
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2)) +
  coord_cartesian(xlim=c(0,1)) +
  scale_y_discrete(labels = label_wrap_gen(width = 22)) 

combined_plot

#ggarrange(combined_plot ,legend_2)
testp <- ggdraw() +
  draw_plot(combined_plot + theme(plot.margin = margin(t = 15, r = 6, b = 6, l = 10, unit = "pt"))
            , x = 0, y = 0, width = 1, height = 0.94) +
  draw_plot(legend_1, x = 0, y = 0.42, width = 1, height = 1) +
  draw_plot(legend_2, x = 0, y = 0.42, width = 1, height = 1)
testp

library(Cairo)

CairoPDF(paste0("./figures/SUPP_switching.pdf"), height = 7, width = 6.5)
print(testp)
dev.off()

svg(paste0("./figures/SUPP_switching.svg"), height = 7, width = 6.5)
print(testp)
dev.off()


# AGO1 PLOT  -----------------------------------------------
#Need to manually run these IsoformSwitchAnalyzeR functions to run the modified function below
source('./code/figures/Functions.R')

  switchAnalyzeRlist = aSwitchList_part2
  gene = "AGO1"
  isoform_id = NULL
  condition1 = "t00"
  condition2 = "t30"
  
  ### Advanced arguments
  IFcutoff = 0.05
  dIFcutoff = 0.05
  alphas = c(0.05, 0.001)
  rescaleTranscripts = TRUE
  plotTopology = FALSE
  reverseMinus = TRUE
  addErrorbars = TRUE
  logYaxis = FALSE
  localTheme = theme_bw(base_size = 6.4)
  additionalArguments  = list(switchPlotTranscript, switchPlotGeneExp, switchPlotIsoExp, switchPlotIsoUsage)
  
  # check conditions
  if (TRUE) {
    ### Identify conditions
    allConditionPairs <-
      unique(switchAnalyzeRlist$isoformFeatures[,c(
        'condition_1', 'condition_2'
      )])
    
    levelsToMatch <- as.vector(t(allConditionPairs))
      
      
      conditionsFound <-
        apply(allConditionPairs, 1, function(x) {
          x[1] == condition1 & x[2] == condition2
        })
      
      if (!any(conditionsFound)) {
        # the the conditions are not found try revers order
        conditionsFound <-
          apply(allConditionPairs, 1, function(x) {
            x[2] == condition1 & x[1] == condition2
          })
        
        # if now found change the order of conditions
        if (any(conditionsFound)) {
          temp <- condition1
          condition1 <- condition2
          condition2 <- temp
        } else {
          stop(
            paste(
              'The wanted comparison: \'',
              condition1,
              ' vs ',
              condition2,
              '\' is not in the data - please revise accordingly'
            )
          )
        }
      }
      
    } else {
      condition1 <- allConditionPairs$condition_1
      condition2 <- allConditionPairs$condition_2
    }
    
    levelsToMatch <-
      levelsToMatch[which(
        levelsToMatch %in% c(condition1, condition2)
      )]
  
  
  ### Parse additional list argument
  if (TRUE) {
    ### Expression plot arguments
    if (TRUE) {
      ## Make list with default values - can be replaced by as.list(args(expressionPlots)) ?
      eArgList <- list(
        switchAnalyzeRlist          = switchAnalyzeRlist,
        gene                        = gene,
        isoform_id                  = isoform_id,
        condition1                  = condition1,
        condition2                  = condition2,
        IFcutoff                    = IFcutoff,
        addErrorbars                = addErrorbars,
        confidenceIntervalErrorbars = TRUE,
        confidenceInterval          = 0.95,
        #alphas                      = c(0.05, 0.001),
        alphas                      = alphas,
        logYaxis                    = logYaxis,
        extendFactor                = 0.05,
        localTheme                  = localTheme
      )
      
      ### Modify default arguments if nessesary
      # subset to only those not already covered by the arguments
      additionalArguments2 <-
        additionalArguments[which(
          names(additionalArguments) %in% c(
            'confidenceIntervalErrorbars',
            'confidenceInterval',
            'alphas',
            'extendFactor',
            'rescaleRoot'
          )
        )]
      if (length(additionalArguments2) != 0) {
        newArgNames <- names(additionalArguments2)
        
        for (i in seq_along(additionalArguments2)) {
          eArgList[[newArgNames[i]]] <- additionalArguments2[[i]]
        }
      }
      
    }
    
  
    
  }

  ### Make subplots - these functions also check the input thoroughly

    # expression plots
    expressionPlots <- expressionAnalysisPlot(
      switchAnalyzeRlist           = eArgList$switchAnalyzeRlist,
      gene                         = eArgList$gene,
      isoform_id                   = eArgList$isoform_id,
      condition1                   = eArgList$condition1,
      condition2                   = eArgList$condition2,
      IFcutoff                     = eArgList$IFcutoff,
      addErrorbars                 = eArgList$addErrorbars,
      confidenceIntervalErrorbars  = eArgList$confidenceIntervalErrorbars,
      confidenceInterval           = eArgList$confidenceInterval,
      alphas                       = eArgList$alphas,
      optimizeForCombinedPlot      = TRUE,
      logYaxis                     = eArgList$logYaxis,
      extendFactor                 = eArgList$extendFactor,
      localTheme                   = eArgList$localTheme
    )
  #   
  
  ### Change plot margin
  marginSize <- 0.1
  expressionPlots2 <-
    lapply(expressionPlots[1:2], function(x) {
      x +
        theme(
          plot.margin = margin(
            t = 6,
            r = 6,
            b = 6,
            l = 6,
            unit = "pt" ),
          legend.position="none"
        )
    })
  
  expressionPlots3 <-
    lapply(expressionPlots[3], function(x) {
      x +
        theme(
          margin(
            t = 6,
            r = 6,
            b = 6,
            l = 6,
            unit = "pt" ),
          legend.position="none"
        )
    })

  label_new <- c(
    "<span style='color:#009E73'><b>PB.1207.309",
    "<span style='color:#009E73'><b>PB.1207.317",
    "<span style='color:#D55E00'><b>PB.1207.336")

#Run this file to bring in the edited AGO1 plotting functions
source("./code/figures/ago_plotting.R")
myPlot <- AGO1_plotting()
  
library(ggtext)
  test <-
    myPlot + 
    theme(plot.margin = margin(
      t = 6,
      r = 6,
      b = 6,
      l = 6,
      unit = "pt"),
      legend.margin = margin(t= -6,0,0,0,unit = "pt"),
      legend.position="bottom", 
      legend.key.size = unit(0.25, 'cm'), 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10), 
      axis.text.y = element_markdown(size = 8), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 8)) +
    scale_y_continuous(
      breaks = c(3,2,1),  # Set the y-axis positions
      labels = label_new  # Apply the custom labels
    ) +
    labs(title = "AGO1 - t00 vs t30") +
    guides(fill = guide_legend(nrow = 1))
  test 
# 
expressionPlots2$gene_expression$data$Analyis <- "DGE"
expressionPlots3$isoform_usage$data$Analyis <- "DTU"
expressionPlots2$isoform_expression$data$Analyis <- "DTE"


AGO1 <- ggdraw() +
  draw_plot(test + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), 
                         plot.title = element_text(hjust = 0.5, size = 8, face = "bold")), x = 0, y = 0.6, width = 1, height = 0.4) +
  draw_plot(expressionPlots2$gene_expression +  
              theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust=1, size = 6.4, colour = "black"), 
                    axis.text.y=element_text(size = 6.4, colour = "black"), 
                    axis.title = element_text(size = 8),
                    strip.text = element_text(size = 8) ),
            x = 0, y = 0, width = 0.3, height = 0.6) +
  draw_plot(expressionPlots3$isoform_usage + 
              theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt"), 
                    axis.text.x=element_text(angle=-45, hjust = 0, vjust=1, size = 6.4, colour = "black"), 
                    axis.text.y=element_text(size = 6.4, colour = "black"), 
                    axis.title = element_text(size = 8),
                    strip.text = element_text(size = 8)) + 
              facet_wrap(~ Analyis, labeller = label_wrap_gen(width = 10)) , x = 0.65, y = 0, width = 0.35, height = 0.6) +
  draw_plot(expressionPlots2$isoform_expression +    
              theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust=1, size = 6.4, colour = "black"), 
                    axis.text.y=element_text(size = 6.4, colour = "black"), 
                    axis.title = element_text(size = 8),
                    strip.text = element_text(size = 8)) +
              facet_wrap(~ Analyis, labeller = label_wrap_gen(width = 10)), x = 0.3, y = 0, width = 0.35, height = 0.6)
AGO1

AGO1 <- AGO1 + theme(plot.margin = margin(t = 0, r = 16, b = 0, l = 6, unit = "pt"))

  

# IR + NMD ----------------------------------------------------------------
#IR transcripts:
AS <- aSwitchList_part2$AlternativeSplicingAnalysis
AS_IR <- AS %>% filter(IR > 0) #12224 transcripts with 1 or more retained introns

AS_IR <- merge(AS_IR, classification[, c("isoform", "associated_gene")],
               by.x = "isoform_id", by.y = "isoform", all.x = T)
CPM <- read.csv("./code/expression/output/transcript_CPM.csv")

CPM$Avgt00 <- rowMeans(subset(CPM, select = c(1:3)))
CPM$Avgt04 <- rowMeans(subset(CPM, select = c(4:9)))
CPM$Avgt30 <- rowMeans(subset(CPM, select = c(10:15)))

CPM$IR <- ifelse(CPM$pb_id %in% AS_IR$isoform_id, "IR", "no_IR")

#Only keep isoforms in the AS from isoswitch object:
CPM_AS <- CPM %>% filter(pb_id %in% aSwitchList_part2$AlternativeSplicingAnalysis$isoform_id)



#PTC transcripts:
PTC_tr <- aSwitchList_part2$isoformFeatures %>% filter(PTC == TRUE) %>% pull(isoform_id) %>% unique()

CPM$PTC <- ifelse(CPM$pb_id %in% PTC_tr, "PTC", "no_PTC")
CPM_PTC <- CPM %>% filter(pb_id %in% aSwitchList_part2$AlternativeSplicingAnalysis$isoform_id)


#IR = 
IR_tr <- CPM %>% filter(IR == "IR") %>% pull(pb_id) %>% unique()

#IR only:
IR_only <- setdiff(IR_tr, PTC_tr) #9606
#NMD only:
PTC_only <- setdiff(PTC_tr, IR_tr)
#Overlap:
overlap <- intersect(IR_tr, PTC_tr)
#Remaining:
all <- unique(union(IR_tr, PTC_tr))
remaining <- setdiff(CPM_AS$pb_id, all)

transcript_groups <- list(
  IR_only = IR_only,
  PTC_only = PTC_only,
  overlap = overlap, 
  remaining = remaining)

# Convert the list into a long data.frame for joining
group_df <- bind_rows(lapply(names(transcript_groups), function(name) {
  data.frame(pb_id = transcript_groups[[name]], group = name)
}))

# Join and label the group
CPM_AS <- CPM_AS %>%
  left_join(group_df, by = "pb_id")


cpm_long_AS <- CPM_AS %>%
  pivot_longer(c("Avgt00", "Avgt04", "Avgt30"), 
               names_to = "Timepoint", 
               values_to = "CPM")

library(rstatix)
library(scales)
cpm_long_AS$log_cpm <- log2(cpm_long_AS$CPM + 1)

stat.test.ir <- cpm_long_AS %>% group_by(Timepoint) %>% wilcox_test(log_cpm ~ group)
stat.test.ir <- stat.test.ir %>% mutate(p.adj.signif = p.adjust(p, method='bonferroni')) #need to do sep, b/c of group_by
stat.test.ir <- stat.test.ir %>% add_xy_position(x = "Timepoint") 

stat.test.ir$p.adj.signif <- sprintf("%.2e", stat.test.ir$p.adj.signif)
stat.test.ir <- stat.test.ir %>% add_significance(p.col = "p.adj")


stat.test.ir$y.position <- log2(stat.test.ir$y.position + 1)

stat.test.ir <- stat.test.ir %>%
  mutate(group1 = recode(group1, "IR_only" = "IR only"), 
         group2 = recode(group2, "IR_only" = "IR only"),
         group1 = recode(group1, "overlap" = "IR & PTC"),
         group2 = recode(group2, "overlap" = "IR & PTC"),
         group1 = recode(group1, "PTC_only" = "PTC only"),
         group2 = recode(group2, "PTC_only" = "PTC only"),
         group1 = recode(group1, "remaining" = "Others"),
         group2 = recode(group2, "remaining" = "Others"))

cpm_long_AS <- cpm_long_AS %>%
  mutate(group = recode(group, "IR_only" = "IR only"), 
         group = recode(group, "IR_only" = "IR only"),
         group = recode(group, "overlap" = "IR & PTC"),
         group = recode(group, "overlap" = "IR & PTC"),
         group = recode(group, "PTC_only" = "PTC only"),
         group = recode(group, "PTC_only" = "PTC only"),
         group = recode(group, "remaining" = "Others"),
         group = recode(group, "remaining" = "Others"))

cpm_long_AS$group <- factor(cpm_long_AS$group, levels = c("IR only", "IR & PTC", "PTC only", "Others"))


bp <- cpm_long_AS %>% 
  ggplot(aes(x=factor(Timepoint), y=log_cpm)) + 
  geom_boxplot(aes(fill = group, colour = group), 
               # outlier.size = 0.0001, 
               outlier.shape = NA, 
               size = 0.2, 
               alpha = 0.35,
               #  size = 0.5, 
               linewidth = 0.45, 
               width = 0.7
  ) +
  ylim(0, 6.4) +
  labs(y = "log2(CPM + 1)") +
  stat_pvalue_manual(stat.test.ir,
                     label = "p.adj.signif",
                     size = 2.5 ,
                     bracket.nudge.y = 1.2,
                     step.group.by = "Timepoint",
                     step.increase = 0.02,
                     tip.length = 0.01, 
                     bracket.size = 0.2) +
  theme_classic(base_size = 8) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        axis.title.x = element_blank(), 
        legend.position = "right", 
        axis.text = element_text(size = 8)) +
  scale_x_discrete(labels = c("Avgt00" = "t00", "Avgt04" = "t04", "Avgt30" = "t30"))
bp
# 
pdf(paste0("./figures/SUPP_IR_NMD_transcript.pdf"), height = 3.5, width = 4)
print(bp)
dev.off()

svg(paste0("./figures/SUPP_IR_NMD_transcript.svg"), height = 3.5, width = 4)
print(bp)
dev.off()



# APA ---------------------------------------------------------------------
#K APA
K_PPAU_t00_v_t30 <- read.csv("./code/APA/output/APA/PPAU_t00_v_t30.csv")
cluster_t30_vs_t00_df <- read.csv("./code/APA/output/APA/Deseq2_cluster_lrmethod_t00_v_t30.csv")
K_PPAU_t00_v_t30 <- merge(K_PPAU_t00_v_t30, cluster_t30_vs_t00_df[, c("X", "log2FoldChange", "padj")],
                by.x = "sequence_cluster", by.y = "X", all.x = T)
K_PPAU_t00_v_t30 <- K_PPAU_t00_v_t30 %>% mutate(Analysis = "Long-read method", 
                                                Comparison = "t00 vs t30")
K_PPAU_t00_v_t30 <- K_PPAU_t00_v_t30[,c(1, 21, 22, 24:27)]
colnames(K_PPAU_t00_v_t30)[1] <- "GeneName"


K_PPAU_t00_v_t04 <- read.csv("./code/APA/output/APA/PPAU_t00_v_t04.csv")
cluster_t04_vs_t00_df <- read.csv("./code/APA/output/APA/Deseq2_cluster_lrmethod_t00_v_t04.csv")
K_PPAU_t00_v_t04 <- merge(K_PPAU_t00_v_t04, cluster_t04_vs_t00_df[, c("X", "log2FoldChange", "padj")],
                          by.x = "sequence_cluster", by.y = "X", all.x = T)
K_PPAU_t00_v_t04 <- K_PPAU_t00_v_t04 %>% mutate(Analysis = "Long-read method", 
                                                Comparison = "t00 vs t04")
K_PPAU_t00_v_t04 <- K_PPAU_t00_v_t04[,c(1, 21, 22, 24:27)]
colnames(K_PPAU_t00_v_t04)[1] <- "GeneName"



K_PPAU_t04_v_t30 <- read.csv("./code/APA/output/APA/PPAU_t04_v_t30.csv")
cluster_t30_vs_t04_df <- read.csv("./code/APA/output/APA/Deseq2_cluster_lrmethod_t04_v_t30.csv")
K_PPAU_t04_v_t30 <- merge(K_PPAU_t04_v_t30, cluster_t30_vs_t04_df[, c("X", "log2FoldChange", "padj")],
                          by.x = "sequence_cluster", by.y = "X", all.x = T)
K_PPAU_t04_v_t30 <- K_PPAU_t04_v_t30 %>% mutate(Analysis = "Long-read method", 
                                                Comparison = "t04 vs t30")
K_PPAU_t04_v_t30 <- K_PPAU_t04_v_t30[,c(1, 21, 22, 24:27)]
colnames(K_PPAU_t04_v_t30)[1] <- "GeneName"




#QAPA
QAPA_PPAU_t00_v_30 <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t00_v_t30.csv")
cluster_t00_vs_t30_df <- read.csv("./code/APA/output/QAPA/Deseq2_cluster_t00_v_t30.csv")
QAPA_PPAU_t00_v_30 <- merge(QAPA_PPAU_t00_v_30, cluster_t00_vs_t30_df[, c("X", "log2FoldChange", "padj")],
                   by.x = "GeneName", by.y = "X", all.x = T)
QAPA_PPAU_t00_v_30 <- QAPA_PPAU_t00_v_30 %>% mutate(Analysis = "QAPA", 
                                                    Comparison = "t00 vs t30")
QAPA_PPAU_t00_v_30 <- QAPA_PPAU_t00_v_30[,c(1, 22:27)]


QAPA_PPAU_t00_v_04 <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t00_v_t04.csv")
cluster_t00_vs_t04_df <- read.csv("./code/APA/output/QAPA/Deseq2_cluster_t00_v_t04.csv")
QAPA_PPAU_t00_v_04 <- merge(QAPA_PPAU_t00_v_04, cluster_t00_vs_t04_df[, c("X", "log2FoldChange", "padj")],
                            by.x = "GeneName", by.y = "X", all.x = T)
QAPA_PPAU_t00_v_04 <- QAPA_PPAU_t00_v_04 %>% mutate(Analysis = "QAPA", 
                                                    Comparison = "t00 vs t04")
QAPA_PPAU_t00_v_04 <- QAPA_PPAU_t00_v_04[,c(1, 22:27)]


QAPA_PPAU_t04_v_30 <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t04_v_t30.csv")
cluster_t04_vs_t30_df <- read.csv("./code/APA/output/QAPA/Deseq2_cluster_t04_v_t30.csv")
QAPA_PPAU_t04_v_30 <- merge(QAPA_PPAU_t04_v_30, cluster_t04_vs_t30_df[, c("X", "log2FoldChange", "padj")],
                            by.x = "GeneName", by.y = "X", all.x = T)
QAPA_PPAU_t04_v_30 <- QAPA_PPAU_t04_v_30 %>% mutate(Analysis = "QAPA", 
                                                    Comparison = "t04 vs t30")
QAPA_PPAU_t04_v_30 <- QAPA_PPAU_t04_v_30[,c(1, 22:27)]




df_list <- list(K_PPAU_t00_v_t04, K_PPAU_t00_v_t30, K_PPAU_t04_v_t30, 
                QAPA_PPAU_t00_v_04, QAPA_PPAU_t00_v_30, QAPA_PPAU_t04_v_30)
df_all <- do.call(rbind, df_list)
df_all <- df_all %>% tidyr::replace_na(list(padj = 1))


library(ggrepel)
library(tidyr)

cluster_plot <- ggplot(df_all, aes(x=dPPAU, y=log2FoldChange)) + 
  geom_point(alpha = 0.5, size = 1, color = ifelse((df_all$log2FoldChange > 1 & df_all$padj <0.05),
                                         "#FF958C",
                                         ifelse((df_all$log2FoldChange < -1 & df_all$padj <0.05), "#883677", "grey"))) + 
#  geom_smooth(method = lm, se = FALSE) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.001, size = 3.2) +
  geom_vline(xintercept = c(-20, 20), linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") + 

  annotate("text", x = 28, y = -6, label = "Lengthening", hjust = 0, size = 2.8) +
  annotate("segment",
           x = 30, xend = 85,
           y = -7, yend = -7,
           arrow = arrow(length = unit(0.2, "cm")),
           colour = "red") +
  
  # Arrow + label for "shorter"
  annotate("text", x = -87, y = -6, label = "Shortening", hjust = 0, size = 2.8) +
  annotate("segment",
           x = -30, xend = -85,
           y = -7, yend = -7,
           arrow = arrow(length = unit(0.2, "cm")),
           colour = "blue") +
  
  theme_bw(base_size = 9) +
  theme(panel.spacing.x = unit(0.56, "cm")) +
  labs(y = "log2FC") +
  facet_grid(factor(Analysis) ~ factor(Comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")) ) 
cluster_plot

pdf(paste0("./figures/SUPP_APA.pdf"), height = 5, width = 6.5)
print(cluster_plot)
dev.off()

svg(paste0("./figures/SUPP_APA.svg"), height = 5, width = 6.5)
print(cluster_plot)
dev.off()



# SFARI GENE ENRICHMENT ---------------------------------------------------
# Function to compute Fisher's exact test for pairwise comparisons
#AS term enrichment:
#Run AS code above (at gene-level):
localConseq4$condition <- paste0(localConseq4$condition_1, "vs", localConseq4$condition_2)
localConseq4_summ <- localConseq4 %>% group_by(condition, AStype, gene_id) %>%
  summarise(Sum = sum(nrDiff)) %>%
  ungroup()

AS <- localConseq4_summ

#All genes (up or down)

#To get all genes:
gene_lists_t00_t04_AS <- AS %>%
  filter(condition == "t00vst04") %>%
  group_by(AStype) %>%
  filter(Sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(AStype, gene_ids) %>%
  tibble::deframe()

gene_lists_t00_t30_AS <- AS %>%
  filter(condition == "t00vst30") %>%
  group_by(AStype) %>%
  filter(Sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(AStype, gene_ids) %>%
  tibble::deframe()

gene_lists_t04_t30_AS <- AS %>%
  filter(condition == "t04vst30") %>%
  group_by(AStype) %>%
  filter(Sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(AStype, gene_ids) %>%
  tibble::deframe()

#Bg genes:
bg.AS <- unique(aSwitchList_part2$isoformFeatures$gene_id) #10,529


SFARI <- read.csv("./data/SFARI-Gene_genes_04-03-2025release_04-15-2025export.csv")


#Testing for enrichment in the each functional consequence:
SFARI.AS <- SFARI %>% filter(gene.symbol %in% bg.AS) %>% pull(gene.symbol) %>% unique()


#Run tests for t00 vs t04:
results_list_t00_t04 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.AS
  setB <- gene_lists_t00_t04_AS[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.AS)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t00_t04_AS[i]))
  results_list_t00_t04[[pair_name]] <- result
}

results_df_t00_t04 <- do.call(rbind, results_list_t00_t04)
results_df_t00_t04 <- data.frame(pair = names(results_list_t00_t04), results_df_t00_t04)
results_df_t00_t04$p_adjusted <- p.adjust(results_df_t00_t04$p_value, method = "bonferroni")

#Run tests for t00 vs t30:
results_list_t00_t30 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.AS
  setB <- gene_lists_t00_t30_AS[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.AS)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t00_t30_AS[i]))
  results_list_t00_t30[[pair_name]] <- result
}

results_df_t00_t30 <- do.call(rbind, results_list_t00_t30)
results_df_t00_t30 <- data.frame(pair = names(results_list_t00_t30), results_df_t00_t30)
results_df_t00_t30$p_adjusted <- p.adjust(results_df_t00_t30$p_value, method = "bonferroni")

#Run tests for t00 vs t04:
results_list_t04_t30 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.AS
  setB <- gene_lists_t04_t30_AS[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.AS)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t04_t30_AS[i]))
  results_list_t04_t30[[pair_name]] <- result
}

results_df_t04_t30 <- do.call(rbind, results_list_t04_t30)
results_df_t04_t30 <- data.frame(pair = names(results_list_t04_t30), results_df_t04_t30)
results_df_t04_t30$p_adjusted <- p.adjust(results_df_t04_t30$p_value, method = "bonferroni")

#Combined data:
results_df_t00_t04$comparison <- "t00 vs t04"
results_df_t00_t30$comparison <- "t00 vs t30"
results_df_t04_t30$comparison <- "t04 vs t30"

results_AS <- rbind(rbind(results_df_t00_t04, results_df_t00_t30), results_df_t04_t30)
results_AS <- results_AS %>%
  mutate(
    consequence = str_remove(pair, "^SFARI vs\\s+")
  )
results_AS

######Func conseq term enrichment:
bg.func <- bg.AS
#SFARI <- read.csv("./data/SFARI-Gene_genes_04-03-2025release_04-15-2025export.csv")


#Testing for enrichment in the each functional consequence:
SFARI.func <- SFARI.AS

#Genes in each conseq:
Conseq_genes_count <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = T,
  analysisOppositeConsequence=F,
  plot=F,
  returnResult=TRUE,
  returnSummary=F)
#Conseq_genes_count <- Conseq_genes_count %>% filter((condition_1 == "t00" & condition_2 == "t30"))
Conseq_genes_count$comparison <- paste0(Conseq_genes_count$condition_1, " vs ", Conseq_genes_count$condition_2)
Conseq_genes_count <- Conseq_genes_count %>% filter(featureCompared == "3_utr_length" |
                                                      featureCompared == "5_utr_length" |
                                                      featureCompared == "domain_length" |
                                                      featureCompared == "domains_identified" |
                                                      featureCompared == "IDR_identified" |
                                                      featureCompared == "IDR_length" |
                                                      featureCompared == "isoform_length" |
                                                      featureCompared == "NMD_status" |
                                                      featureCompared == "ORF_length" 
                                                     # featureCompared == "tss" |
                                                    #  featureCompared == "tts"
                                                      )

Conseq_genes_count <- Conseq_genes_count %>% group_by(comparison, gene_id) %>%
                                       mutate(nrDiff = ifelse((switchConsequence == "3UTR is longer"|
                                       switchConsequence == "5UTR is longer"|
                                       switchConsequence == "Domain gain"|
                                       switchConsequence == "Domain length gain"|
                                       switchConsequence == "IDR gain"|
                                       switchConsequence == "IDR length gain"|
                                       switchConsequence == "Length gain"|
                                       switchConsequence == "NMD sensitive"| ##This is the opposite.
                                       switchConsequence == "ORF is longer"
                                      # switchConsequence == "Tss more upstream"|
                                      # switchConsequence == "Tts more downstream"
                                        ), 1, -1)) %>% ungroup()
info <- Conseq_genes_count %>% dplyr::select(switchConsequence, featureCompared) %>% unique()

Conseq_genes_summ <- Conseq_genes_count %>% group_by(comparison, switchConsequence, gene_id) %>%
  dplyr::summarize(sum = sum(nrDiff)) %>%
  ungroup() 
Conseq_genes_summ <- merge(Conseq_genes_summ, info[,c('switchConsequence', 'featureCompared')])


#To get all genes:
gene_lists_t00_t04 <- Conseq_genes_summ %>%
  filter(comparison == "t00 vs t04") %>%
  group_by(featureCompared) %>%
  filter(sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(featureCompared, gene_ids) %>%
  tibble::deframe()
  
gene_lists_t00_t30 <- Conseq_genes_summ %>%
  filter(comparison == "t00 vs t30") %>%
  group_by(featureCompared) %>%
  filter(sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(featureCompared, gene_ids) %>%
  tibble::deframe()

gene_lists_t04_t30 <- Conseq_genes_summ %>%
  filter(comparison == "t04 vs t30") %>%
  group_by(featureCompared) %>%
  filter(sum != 0) %>%
  summarise(gene_ids = list(unique(gene_id)), .groups = "drop") %>%
  dplyr::select(featureCompared, gene_ids) %>%
  tibble::deframe()

#Run tests for t00 vs t04:
results_list_t00_t04 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.func
  setB <- gene_lists_t00_t04[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.func)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t00_t04[i]))
  results_list_t00_t04[[pair_name]] <- result
}

results_df_t00_t04 <- do.call(rbind, results_list_t00_t04)
results_df_t00_t04 <- data.frame(pair = names(results_list_t00_t04), results_df_t00_t04)
results_df_t00_t04$p_adjusted <- p.adjust(results_df_t00_t04$p_value, method = "bonferroni")

#Run tests for t00 vs t30:
results_list_t00_t30 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.func
  setB <- gene_lists_t00_t30[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.func)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t00_t30[i]))
  results_list_t00_t30[[pair_name]] <- result
}

results_df_t00_t30 <- do.call(rbind, results_list_t00_t30)
results_df_t00_t30 <- data.frame(pair = names(results_list_t00_t30), results_df_t00_t30)
results_df_t00_t30$p_adjusted <- p.adjust(results_df_t00_t30$p_value, method = "bonferroni")

#Run tests for t00 vs t04:
results_list_t04_t30 <- list()

i = 1
for (i in 1:9) {
  setA <- SFARI.func
  setB <- gene_lists_t04_t30[[i]]
  result <- compute_fisher_enrich(setA, setB, bg.func)
  
  # Store the results in the list, and create the pair name
  pair_name <- paste0("SFARI vs ", names(gene_lists_t04_t30[i]))
  results_list_t04_t30[[pair_name]] <- result
}

results_df_t04_t30 <- do.call(rbind, results_list_t04_t30)
results_df_t04_t30 <- data.frame(pair = names(results_list_t04_t30), results_df_t04_t30)
results_df_t04_t30$p_adjusted <- p.adjust(results_df_t04_t30$p_value, method = "bonferroni")

#Combined data:
results_df_t00_t04$comparison <- "t00 vs t04"
results_df_t00_t30$comparison <- "t00 vs t30"
results_df_t04_t30$comparison <- "t04 vs t30"

results_func_conseq <- rbind(rbind(results_df_t00_t04, results_df_t00_t30), results_df_t04_t30)
results_func_conseq <- results_func_conseq %>%
  mutate(
    consequence = str_remove(pair, "^SFARI vs\\s+")
  )

results_AS$conseq_type <- "Alternative Splicing"
results_func_conseq$conseq_type <- "Functional Consequence"

###Combined results + prep:
combined_enrich <- rbind(results_AS, results_func_conseq)

combined_enrich$consequence <- str_replace(combined_enrich$consequence, "MES",
                                           "Multiple exon splicing")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "3_utr_length",
                                     "3' UTR length")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "5_utr_length",
                                     "5' UTR length")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "domain_length",
                                     "Domain length")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "domains_identified",
                                     "Domain gain/loss")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "IDR_identified",
                                     "IDR gain/loss")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "IDR_length",
                                               "IDR length")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "isoform_length",
                                               "Isoform length")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "NMD_status",
                                               "NMD status")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "ORF_length",
                                               "ORF length")

combined_enrich$consequence <- str_replace(combined_enrich$consequence, "A3",
                                               "Alternative 3' splice site")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "A5",
                                               "Alternative 5' splice site")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "ES",
                                               "Exon splicing")
combined_enrich$consequence <- str_replace(combined_enrich$consequence, "IR",
                                               "Intron splicing")

combined_enrich$consequence <- str_replace(combined_enrich$consequence, "MIC",
                                               "Microexon splicing")



combined_enrich <- combined_enrich %>%
  mutate(consequence = fct_relevel(consequence, 
                                   "Alternative 3' splice site", 
                                   "Alternative 5' splice site", 
                                   "Exon splicing",
                                   "Microexon splicing", 
                                   "Intron splicing", 
                                   "Multiple exon splicing",
                               "IDR length", 
                               "IDR gain/loss", 
                               "Domain length", 
                               "Domain gain/loss", 
                               "NMD status", 
                               "ORF length", 
                               "Isoform length",
                               "3' UTR length", 
                               "5' UTR length"))

combined_enrich$sig <- ifelse(combined_enrich$p_adjusted <= 0.05, TRUE, FALSE)
combined_enrich$comparison <- factor(combined_enrich$comparison, levels = c("t00 vs t04", "t04 vs t30", "t00 vs t30"))


##FIg:
combined_enrich_filter <- combined_enrich %>% filter(consequence != "ATTS" &
                                                       consequence != "ATSS" &
                                                       consequence != "MEE" &
                                                       consequence != "Isoform length" &
                                                       consequence != "Domain length" &
                                                       consequence != "IDR length" &
                                                       consequence != "Alternative 3' splice site" &
                                                       consequence != "Alternative 5' splice site" &
                                                       consequence != "Multiple exon splicing" &
                                                       consequence != "ORF length")

desired_order <- c("Exon splicing",
                   "Microexon splicing",
                   "Intron splicing",
                   "NMD status",
                   "Domain gain/loss", 
                   "IDR gain/loss", 
                   "5' UTR length",
                   "3' UTR length")
combined_enrich_filter <- combined_enrich_filter %>%
  mutate(consequence = factor(consequence, levels = rev(desired_order)))
# 

##Plot:
ASD_bar <- combined_enrich_filter %>%
  ggplot(aes(x = odds_ratio.odds.ratio, 
             y = consequence, 
             fill = overlap)) +
  geom_col(width = 0.6) + 
  facet_wrap(~ factor(comparison, c("t00 vs t04", "t04 vs t30", "t00 vs t30")),
             scales = "free_x") +
  scale_fill_gradient(name = "Gene Count", low = "lightblue", high = "darkblue") +
  scale_x_continuous(limits = c(0, 2.85),
                     expand = expansion(mult = c(0, .1))) +
  labs(title = "SFARI Gene Enrichment",
       x = "Odds Ratio",
       y = NULL) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "black") +
  geom_text(data = subset(combined_enrich_filter, sig == TRUE),
            aes(x = odds_ratio.odds.ratio, 
                y = consequence,
                label = "*"),
            hjust = -0.4,     
            vjust = 0.75,
            colour = "black",
            size = 3.5) +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, 
                              size = 8, face = "bold"),
    strip.text = element_text(size = 8),
    axis.text = element_text(colour = "black", size = 6.4),
    axis.title.y = element_blank(), 
    legend.position = "none"
    ) 
ASD_bar

legend_ASD_bar <- ggpubr::get_legend(ASD_bar +
                                 theme( legend.position = "top", 
                                        legend.box = "horizontal",
                                    
                                        legend.text = element_text(size = 7),
                                        legend.key.size = unit(0.25, 'cm'), 
                                      
                                        legend.margin = margin(
                                          t = 0,
                                          r = 0,
                                          b = 0,
                                          l = 0,
                                          unit = "cm"
                                        ) )+
                                          guides(fill = guide_colourbar(title.vjust = 0.8, 
                                                                        barwidth = 3.5))
                                        )

combined_ASD_bar <- ggdraw() +
  draw_plot(ASD_bar + theme(plot.title = element_text(vjust = 16), 
                            axis.title.x = element_text(vjust = -6))
            , x = 0, y = 0, width = 1, height = 0.934) +
  draw_plot(legend_ASD_bar, x = 0, y = 0.447, width = 1, height = 1)
combined_ASD_bar


#

# GO TERMS - SWITCHING GENES ----------------------------------------------
library(gprofiler2)
run_GO_enrichment <- function(gene_list, bg, searches) {
  go_result <- gost(gene_list, 
                    organism = "hsapiens", 
                    significant = T,
                    exclude_iea = F, 
                    user_threshold = 0.05,
                    correction_method = "fdr", 
                    custom_bg = bg, 
                    sources = searches)
  return(go_result)
}

#Run AS code above (at gene-level):
localConseq4$condition <- paste0(localConseq4$condition_1, "vs", localConseq4$condition_2)
localConseq4_summ <- localConseq4 %>% group_by(condition, AStype, gene_id) %>%
  summarise(Sum = sum(nrDiff)) %>%
  ungroup()

###T00 vs T30
AS_t00_t30 <- localConseq4_summ %>% filter(condition == "t00vst30")

#All genes (up or down)
genes_ATTS <- AS_t00_t30 %>% filter(AStype == "ATTS" & Sum != 0) %>% pull(gene_id)
genes_MIC <- AS_t00_t30 %>% filter(AStype == "MIC" & Sum != 0) %>% pull(gene_id)
genes_A3 <- AS_t00_t30 %>% filter(AStype == "A3" & Sum != 0) %>% pull(gene_id)
genes_A5 <- AS_t00_t30 %>% filter(AStype == "A5" & Sum != 0) %>% pull(gene_id)
genes_ATSS <- AS_t00_t30 %>% filter(AStype == "ATSS" & Sum != 0) %>% pull(gene_id)
genes_ES <- AS_t00_t30 %>% filter(AStype == "ES" & Sum != 0) %>% pull(gene_id)
genes_MES <- AS_t00_t30 %>% filter(AStype == "MES" & Sum != 0) %>% pull(gene_id)
genes_IR <- AS_t00_t30 %>% filter(AStype == "IR" & Sum != 0) %>% pull(gene_id)

#Bg genes:
all_AS_genes <- unique(aSwitchList_part2$isoformFeatures$gene_id) #10529


#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
gene_lists <- list(genes_ATTS = genes_ATTS,
                   genes_MIC = genes_MIC,
                   genes_A3 = genes_A3,
                   genes_A5 = genes_A5,
                   genes_ATSS = genes_ATSS,
                   genes_ES = genes_ES,
                   genes_MES = genes_MES,
                   genes_IR = genes_IR)

# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = all_AS_genes,
                                                                  searches = searches)
}
# 

go_results_IR <- enrichment_results$genes_IR$result
go_results_IR$query <- "IR GO Terms - t00 vs t30"
# RP_IR <- gene_lists[["genes_IR"]][grepl("^RP", gene_lists[["genes_IR"]])]
go_results_IR <- go_results_IR %>%
  group_by(source) %>%
  arrange(p_value) %>%
  dplyr::slice(1:5) %>%
  ungroup()
go_results_IR$term_name_wrapped <- str_wrap(go_results_IR$term_name, width = 20)
# Create a dot plot for each gene list and GO term
IR_go <- ggplot(go_results_IR, aes(x = reorder(term_name_wrapped, -p_value), y = -log10(p_value), fill = intersection_size)) +
  geom_col(width = c(0.6*1, 0.6*1.1, 0.6*1.1, 0.6*1.1, 0.6*1.1, 0.6*1.1, 0.6*1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "IR GO Terms - t00 vs t30",
       y = expression(paste(-log[10], "(adjusted ", italic("P"), ")")) ) +
  theme_bw(base_size = 8) +
  coord_flip() +
  scale_fill_gradient(name = "Gene Count", low = "lightblue", high = "darkblue") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, colour = "black", face = "bold"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 6.4, colour = "black"),
        axis.title = element_text(size = 8, colour = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 6.4, colour = "black"),
        #  plot.margin = margin(t = 6, b = 0, l = 16, r = 6, unit = "pt"),
        legend.margin = margin(t = -5, b = 0, l = -25, r = 0, unit = "pt"),
        legend.key.size = unit(0.25, 'cm')
  ) +
  facet_wrap(~source,
             scales = "free_y",
             dir = "v",
             strip.position = "right",
             ncol = 1) +  # Facet by source if applicable
  force_panelsizes(row = c(1.25, 4.5, 1.25))
IR_go




# EXPORT FIG 4 ------------------------------------------------------------
library(Cairo)

Figure_4 <- ggdraw() +
  draw_plot(bars_by_tr, x = 0, y = 0.78, width = 0.5, height = 0.22) +
  draw_plot(final_plot, x = 0.5, y = 0.78, width = 0.5, height = 0.22) +
  draw_plot(testp, x = 0, y = 0.35, width = 0.5, height = 0.43) +
  draw_plot(combined_ASD_bar, x = 0.5, y = 0.38, width = 0.5, height = 0.3815) +
  draw_plot(IR_go, x = 0.01, y = 0, width = 0.35, height = 0.35) +
  draw_plot(AGO1, x = 0.35, y = 0, width = 0.65, height = 0.35) +
  
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), 
                  size = 12, 
                  x = c(0, 0.5, 0, 0.5, 0, 0.355),
                  y = c(1, 1, 0.78, 0.78, 0.35, 0.35))
Figure_4
# 

CairoPDF(paste0("./figures/Fig_4.pdf"), height = 8, width = 6.5)
print(Figure_4)
dev.off()


svg(paste0("./figures/Fig_4.svg"), height = 8, width = 6.5)
print(Figure_4)
dev.off()





# SUPP. FIGS & TABLES ------------------------------------------------------


# GO TERMS - SWITCHING GENES ----------------------------------------------
library(gprofiler2)
run_GO_enrichment <- function(gene_list, bg, searches) {
  # Convert gene symbols to Entrez IDs (if necessary)
  #  gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  go_result <- gost(gene_list, 
                    organism = "hsapiens", 
                    significant = T,
                    exclude_iea = F, 
                    user_threshold = 0.05,
                    correction_method = "fdr", 
                    custom_bg = bg, 
                    sources = searches)
  return(go_result)
}


#Run AS code above (at gene-level):
localConseq4$condition <- paste0(localConseq4$condition_1, "vs", localConseq4$condition_2)
localConseq4_summ <- localConseq4 %>% group_by(condition, AStype, gene_id) %>%
  summarise(Sum = sum(nrDiff)) %>%
  ungroup()

###T00 vs T30
AS_t00_t30 <- localConseq4_summ %>% filter(condition == "t00vst30")

#All genes (up or down)
genes_ATTS <- AS_t00_t30 %>% filter(AStype == "ATTS" & Sum != 0) %>% pull(gene_id)
genes_MIC <- AS_t00_t30 %>% filter(AStype == "MIC" & Sum != 0) %>% pull(gene_id)
genes_A3 <- AS_t00_t30 %>% filter(AStype == "A3" & Sum != 0) %>% pull(gene_id)
genes_A5 <- AS_t00_t30 %>% filter(AStype == "A5" & Sum != 0) %>% pull(gene_id)
genes_ATSS <- AS_t00_t30 %>% filter(AStype == "ATSS" & Sum != 0) %>% pull(gene_id)
genes_ES <- AS_t00_t30 %>% filter(AStype == "ES" & Sum != 0) %>% pull(gene_id)
genes_MES <- AS_t00_t30 %>% filter(AStype == "MES" & Sum != 0) %>% pull(gene_id)
genes_IR <- AS_t00_t30 %>% filter(AStype == "IR" & Sum != 0) %>% pull(gene_id)

#Bg genes:
all_AS_genes <- unique(aSwitchList_part2$isoformFeatures$gene_id) #10529


#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
gene_lists <- list(genes_ATTS = genes_ATTS,
                   genes_MIC = genes_MIC,
                   genes_A3 = genes_A3,
                   genes_A5 = genes_A5,
                   genes_ATSS = genes_ATSS,
                   genes_ES = genes_ES,
                   genes_MES = genes_MES,
                   genes_IR = genes_IR)

# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = all_AS_genes,
                                                                  searches = searches)
}
# 

###T00 vs T04
AS_t00_t04 <- localConseq4_summ %>% filter(condition == "t00vst04")

#All genes (up or down)
genes_ATTS <- AS_t00_t04 %>% filter(AStype == "ATTS" & Sum != 0) %>% pull(gene_id)
genes_MIC <- AS_t00_t04 %>% filter(AStype == "MIC" & Sum != 0) %>% pull(gene_id)
genes_A3 <- AS_t00_t04 %>% filter(AStype == "A3" & Sum != 0) %>% pull(gene_id)
genes_A5 <- AS_t00_t04 %>% filter(AStype == "A5" & Sum != 0) %>% pull(gene_id)
genes_ATSS <- AS_t00_t04 %>% filter(AStype == "ATSS" & Sum != 0) %>% pull(gene_id)
genes_ES <- AS_t00_t04 %>% filter(AStype == "ES" & Sum != 0) %>% pull(gene_id)
genes_MES <- AS_t00_t04 %>% filter(AStype == "MES" & Sum != 0) %>% pull(gene_id)
genes_IR <- AS_t00_t04 %>% filter(AStype == "IR" & Sum != 0) %>% pull(gene_id)

#Bg genes:
all_AS_genes 

#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
gene_lists <- list(genes_ATTS = genes_ATTS,
                   genes_MIC = genes_MIC,
                   genes_A3 = genes_A3,
                   genes_A5 = genes_A5,
                   genes_ATSS = genes_ATSS,
                   genes_ES = genes_ES,
                   genes_MES = genes_MES,
                   genes_IR = genes_IR)

# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = all_AS_genes,
                                                                  searches = searches)
}


###T04 vs T30
AS_t04_t30 <- localConseq4_summ %>% filter(condition == "t04vst30")

#All genes (up or down)
genes_ATTS <- AS_t04_t30 %>% filter(AStype == "ATTS" & Sum != 0) %>% pull(gene_id)
genes_MIC <- AS_t04_t30 %>% filter(AStype == "MIC" & Sum != 0) %>% pull(gene_id)
genes_A3 <- AS_t04_t30 %>% filter(AStype == "A3" & Sum != 0) %>% pull(gene_id)
genes_A5 <- AS_t04_t30 %>% filter(AStype == "A5" & Sum != 0) %>% pull(gene_id)
genes_ATSS <- AS_t04_t30 %>% filter(AStype == "ATSS" & Sum != 0) %>% pull(gene_id)
genes_ES <- AS_t04_t30 %>% filter(AStype == "ES" & Sum != 0) %>% pull(gene_id)
genes_MES <- AS_t04_t30 %>% filter(AStype == "MES" & Sum != 0) %>% pull(gene_id)
genes_IR <- AS_t04_t30 %>% filter(AStype == "IR" & Sum != 0) %>% pull(gene_id)

#Bg genes:
all_AS_genes 

#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
gene_lists <- list(genes_ATTS = genes_ATTS,
                   genes_MIC = genes_MIC,
                   genes_A3 = genes_A3,
                   genes_A5 = genes_A5,
                   genes_ATSS = genes_ATSS,
                   genes_ES = genes_ES,
                   genes_MES = genes_MES,
                   genes_IR = genes_IR)

# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = all_AS_genes,
                                                                  searches = searches)
}



# GO TERMS - FUNC CONSEQ --------------------------------------------------
###Do for func. conseq:
Conseq_genes <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = T,
  analysisOppositeConsequence=F,
  plot=F,
  returnResult=TRUE,
  returnSummary=T)


##T00 vs T30:
Conseq_genes_count <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = T,
  analysisOppositeConsequence=F,
  plot=F,
  returnResult=TRUE,
  returnSummary=F)
Conseq_genes_count <- Conseq_genes_count %>% filter((condition_1 == "t00" & condition_2 == "t30"))

Conseq_genes_count$nrDiff <- ifelse((Conseq_genes_count$switchConsequence == "3UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "5UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "Domain gain"|
                                       Conseq_genes_count$switchConsequence == "Domain length gain"|
                                       Conseq_genes_count$switchConsequence == "IDR gain"|
                                       Conseq_genes_count$switchConsequence == "IDR length gain"|
                                       Conseq_genes_count$switchConsequence == "Length gain"|
                                       Conseq_genes_count$switchConsequence == "NMD sensitive"| ##This is the opposite.
                                       Conseq_genes_count$switchConsequence == "ORF is longer"|
                                       Conseq_genes_count$switchConsequence == "Tss more upstream"|
                                       Conseq_genes_count$switchConsequence == "Tts more downstream"), 1, -1)

Conseq_genes_summ <- Conseq_genes_count %>% 
  group_by(switchConsequence, gene_id) %>%
  summarize(sum = sum(nrDiff)) %>%
  ungroup()
summary_genes <- Conseq_genes_summ %>% 
  group_by(switchConsequence) %>%
  filter(sum != 0)
# 

#I can just take the genes that match the ones above (b/c all others are 0)
genes_3UTR <- Conseq_genes_summ %>% filter(switchConsequence == "3UTR is longer"  |
                                             switchConsequence == "3UTR is shorter") %>% pull(gene_id)
genes_5UTR <- Conseq_genes_summ %>% filter(switchConsequence == "5UTR is longer" |
                                             switchConsequence == "5UTR is shorter") %>% pull(gene_id)
genes_domain_number <- Conseq_genes_summ %>% filter(switchConsequence == "Domain gain" |
                                                      switchConsequence == "Domain loss") %>% pull(gene_id)
genes_domain_length <- Conseq_genes_summ %>% filter(switchConsequence == "Domain length gain" |
                                                      switchConsequence == "Domain length loss") %>% pull(gene_id)
genes_IDR_number <- Conseq_genes_summ %>% filter(switchConsequence == "IDR gain" |
                                                   switchConsequence == "IDR loss") %>% pull(gene_id)
genes_IDR_length <- Conseq_genes_summ %>% filter(switchConsequence == "IDR length gain" |
                                                   switchConsequence == "IDR length loss") %>% pull(gene_id)
genes_length <- Conseq_genes_summ %>% filter(switchConsequence == "Length gain" |
                                               switchConsequence == "Length loss") %>% pull(gene_id)
genes_NMD <- Conseq_genes_summ %>% filter(switchConsequence == "NMD sensitive" |
                                            switchConsequence == "NMD insensitive") %>% pull(gene_id)
genes_ORF_length <- Conseq_genes_summ %>% filter(switchConsequence == "ORF is longer" |
                                                   switchConsequence == "ORF is shorter") %>% pull(gene_id)
genes_TSS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tss more upstream" |
                                                   switchConsequence == "Tss more downstream") %>% pull(gene_id)
genes_TTS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tts more downstream" |
                                                   switchConsequence == "Tts more upstream") %>% pull(gene_id)


# Your gene lists in a named list
gene_lists <- list(
  genes_3UTR = genes_3UTR,
  genes_5UTR = genes_5UTR,
  genes_domain_number = genes_domain_number, 
  genes_domain_length = genes_domain_length, 
  genes_IDR_number = genes_IDR_number, 
  genes_IDR_length = genes_IDR_length, 
  genes_length = genes_length, 
  genes_NMD = genes_NMD, 
  genes_ORF_length = genes_ORF_length, 
  genes_TSS_change = genes_TSS_change, 
  genes_TTS_change = genes_TTS_change)


#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
bg_genes <- unique(aSwitchList_part2$isoformFeatures$gene_id)
# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = bg_genes,
                                                                  searches = searches)
}
#go_results_IDR_number <- enrichment_results$genes_IDR_number$result
t00_t30 <- enrichment_results



##T00 vs T04:
Conseq_genes_count <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = T,
  analysisOppositeConsequence=F,
  plot=F,
  returnResult=TRUE,
  returnSummary=F)
Conseq_genes_count <- Conseq_genes_count %>% filter((condition_1 == "t00" & condition_2 == "t04"))

Conseq_genes_count$nrDiff <- ifelse((Conseq_genes_count$switchConsequence == "3UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "5UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "Domain gain"|
                                       Conseq_genes_count$switchConsequence == "Domain length gain"|
                                       Conseq_genes_count$switchConsequence == "IDR gain"|
                                       Conseq_genes_count$switchConsequence == "IDR length gain"|
                                       Conseq_genes_count$switchConsequence == "Length gain"|
                                       Conseq_genes_count$switchConsequence == "NMD sensitive"| ##This is the opposite.
                                       Conseq_genes_count$switchConsequence == "ORF is longer"|
                                       Conseq_genes_count$switchConsequence == "Tss more upstream"|
                                       Conseq_genes_count$switchConsequence == "Tts more downstream"), 1, -1)

Conseq_genes_summ <- Conseq_genes_count %>% 
  group_by(switchConsequence, gene_id) %>%
  summarize(sum = sum(nrDiff)) %>%
  ungroup()
summary_genes <- Conseq_genes_summ %>% 
  group_by(switchConsequence) %>%
  filter(sum != 0)

#I can just take the genes that match the ones above (b/c all others are 0)
genes_3UTR <- Conseq_genes_summ %>% filter(switchConsequence == "3UTR is longer"  |
                                             switchConsequence == "3UTR is shorter") %>% pull(gene_id)
genes_5UTR <- Conseq_genes_summ %>% filter(switchConsequence == "5UTR is longer" |
                                             switchConsequence == "5UTR is shorter") %>% pull(gene_id)
genes_domain_number <- Conseq_genes_summ %>% filter(switchConsequence == "Domain gain" |
                                                      switchConsequence == "Domain loss") %>% pull(gene_id)
genes_domain_length <- Conseq_genes_summ %>% filter(switchConsequence == "Domain length gain" |
                                                      switchConsequence == "Domain length loss") %>% pull(gene_id)
genes_IDR_number <- Conseq_genes_summ %>% filter(switchConsequence == "IDR gain" |
                                                   switchConsequence == "IDR loss") %>% pull(gene_id)
genes_IDR_length <- Conseq_genes_summ %>% filter(switchConsequence == "IDR length gain" |
                                                   switchConsequence == "IDR length loss") %>% pull(gene_id)
genes_length <- Conseq_genes_summ %>% filter(switchConsequence == "Length gain" |
                                               switchConsequence == "Length loss") %>% pull(gene_id)
genes_NMD <- Conseq_genes_summ %>% filter(switchConsequence == "NMD sensitive" |
                                            switchConsequence == "NMD insensitive") %>% pull(gene_id)
genes_ORF_length <- Conseq_genes_summ %>% filter(switchConsequence == "ORF is longer" |
                                                   switchConsequence == "ORF is shorter") %>% pull(gene_id)
genes_TSS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tss more upstream" |
                                                   switchConsequence == "Tss more downstream") %>% pull(gene_id)
genes_TTS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tts more downstream" |
                                                   switchConsequence == "Tts more upstream") %>% pull(gene_id)


# Your gene lists in a named list
gene_lists <- list(
  genes_3UTR = genes_3UTR,
  genes_5UTR = genes_5UTR,
  genes_domain_number = genes_domain_number, 
  genes_domain_length = genes_domain_length, 
  genes_IDR_number = genes_IDR_number, 
  genes_IDR_length = genes_IDR_length, 
  genes_length = genes_length, 
  genes_NMD = genes_NMD, 
  genes_ORF_length = genes_ORF_length, 
  genes_TSS_change = genes_TSS_change, 
  genes_TTS_change = genes_TTS_change)


#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
bg_genes
# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = bg_genes,
                                                                  searches = searches)
}
t00_t04 <- enrichment_results



##T04 vs T30:
Conseq_genes_count <- extractConsequenceEnrichment(
  aSwitchList_part2,
  alpha=0.05,
  dIFcutoff = 0.05,
  countGenes = T,
  analysisOppositeConsequence=F,
  plot=F,
  returnResult=TRUE,
  returnSummary=F)
Conseq_genes_count <- Conseq_genes_count %>% filter((condition_1 == "t04" & condition_2 == "t30"))

Conseq_genes_count$nrDiff <- ifelse((Conseq_genes_count$switchConsequence == "3UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "5UTR is longer"|
                                       Conseq_genes_count$switchConsequence == "Domain gain"|
                                       Conseq_genes_count$switchConsequence == "Domain length gain"|
                                       Conseq_genes_count$switchConsequence == "IDR gain"|
                                       Conseq_genes_count$switchConsequence == "IDR length gain"|
                                       Conseq_genes_count$switchConsequence == "Length gain"|
                                       Conseq_genes_count$switchConsequence == "NMD sensitive"| ##This is the opposite.
                                       Conseq_genes_count$switchConsequence == "ORF is longer"|
                                       Conseq_genes_count$switchConsequence == "Tss more upstream"|
                                       Conseq_genes_count$switchConsequence == "Tts more downstream"), 1, -1)

Conseq_genes_summ <- Conseq_genes_count %>% 
  group_by(switchConsequence, gene_id) %>%
  summarize(sum = sum(nrDiff)) %>%
  ungroup()
summary_genes <- Conseq_genes_summ %>% 
  group_by(switchConsequence) %>%
  filter(sum != 0)


#I can just take the genes that match the ones above (b/c all others are 0)
genes_3UTR <- Conseq_genes_summ %>% filter(switchConsequence == "3UTR is longer"  |
                                             switchConsequence == "3UTR is shorter") %>% pull(gene_id)
genes_5UTR <- Conseq_genes_summ %>% filter(switchConsequence == "5UTR is longer" |
                                             switchConsequence == "5UTR is shorter") %>% pull(gene_id)
genes_domain_number <- Conseq_genes_summ %>% filter(switchConsequence == "Domain gain" |
                                                      switchConsequence == "Domain loss") %>% pull(gene_id)
genes_domain_length <- Conseq_genes_summ %>% filter(switchConsequence == "Domain length gain" |
                                                      switchConsequence == "Domain length loss") %>% pull(gene_id)
genes_IDR_number <- Conseq_genes_summ %>% filter(switchConsequence == "IDR gain" |
                                                   switchConsequence == "IDR loss") %>% pull(gene_id)
genes_IDR_length <- Conseq_genes_summ %>% filter(switchConsequence == "IDR length gain" |
                                                   switchConsequence == "IDR length loss") %>% pull(gene_id)
genes_length <- Conseq_genes_summ %>% filter(switchConsequence == "Length gain" |
                                               switchConsequence == "Length loss") %>% pull(gene_id)
genes_NMD <- Conseq_genes_summ %>% filter(switchConsequence == "NMD sensitive" |
                                            switchConsequence == "NMD insensitive") %>% pull(gene_id)
genes_ORF_length <- Conseq_genes_summ %>% filter(switchConsequence == "ORF is longer" |
                                                   switchConsequence == "ORF is shorter") %>% pull(gene_id)
genes_TSS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tss more upstream" |
                                                   switchConsequence == "Tss more downstream") %>% pull(gene_id)
genes_TTS_change <- Conseq_genes_summ %>% filter(switchConsequence == "Tts more downstream" |
                                                   switchConsequence == "Tts more upstream") %>% pull(gene_id)


# Your gene lists in a named list
gene_lists <- list(
  genes_3UTR = genes_3UTR,
  genes_5UTR = genes_5UTR,
  genes_domain_number = genes_domain_number, 
  genes_domain_length = genes_domain_length, 
  genes_IDR_number = genes_IDR_number, 
  genes_IDR_length = genes_IDR_length, 
  genes_length = genes_length, 
  genes_NMD = genes_NMD, 
  genes_ORF_length = genes_ORF_length, 
  genes_TSS_change = genes_TSS_change, 
  genes_TTS_change = genes_TTS_change)


#GO terms for all AS results
searches <- c("GO:BP", "GO:MF", "GO:CC", "MIRNA")
bg_genes 
# Initialize an empty list to store results
enrichment_results <- list()
i = 1
# Loop through each gene list and perform GO enrichment
for (i in seq_along(gene_lists)) {
  enrichment_results[[names(gene_lists)[i]]] <- run_GO_enrichment(gene_lists[[i]],
                                                                  bg = bg_genes,
                                                                  searches = searches)
}
t04_t30 <- enrichment_results


