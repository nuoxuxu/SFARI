# ATSS-AS COORDINATION

library(readr)       
library(dplyr)       
library(stringr)     
library(tibble)     
library(tidyr)       
library(forcats)     
library(ggplot2)     
library(cowplot)     
library(patchwork)   
library(purrr)
library(arrow)
library(emmeans)    
library(scales)
library(ggtext)
library(GenomicRanges)

library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(readr)


# --- Load Raw Data ---
raw_data <- read_csv("./code/input/TSS_counts.csv", show_col_types = FALSE)




#Function to filter TSSs
filter_TSS_sites_by_replicates_new <- function(data_frame, min_count = 5, min_replicates = 3, min_count_all_reps = 1) {
  
  # Step 1: Filter TSS_ids that meet min_count_all_reps in *every* sample (ignoring type)
  if (!is.null(min_count_all_reps)) {
    TSS_ids_to_keep <- data_frame %>%
      group_by(event_id, TSS_id, timepoint) %>%
      summarise(total_sum = sum(total, na.rm = TRUE), .groups = "drop") %>%
      group_by(event_id, TSS_id) %>%
      summarise(all_timepoints_ok = all(total_sum >= min_count_all_reps), .groups = "drop") %>%
      filter(all_timepoints_ok) %>%
      dplyr::select(event_id, TSS_id)
    
    data_frame <- data_frame %>%
      semi_join(TSS_ids_to_keep, by = c("event_id", "TSS_id"))
    
    cat("Retained", nrow(TSS_ids_to_keep), "pA sites passing min_count_all_reps =", min_count_all_reps, "\n")
  }
  # Step 2: Filter rows that meet the read count threshold
  filtered <- data_frame %>%
    filter(total >= min_count)
  
  # Step 3: For each event_id and TSS_id, count how many unique timepoints meet the threshold
  sites_to_keep <- filtered %>%
    group_by(event_id, TSS_id) %>%
    summarise(
      n_replicates = n_distinct(timepoint),
      .groups = 'drop'
    ) %>%
    filter(n_replicates >= min_replicates)
  
  cat("Filtered to", nrow(sites_to_keep), "TSS sites across", 
      n_distinct(sites_to_keep$event_id), "event IDs meeting thresholds.\n")
  
  # Step 4: Use semi_join to keep only the matching rows from the original data
  filtered_data <- data_frame %>%
    semi_join(sites_to_keep, by = c("event_id", "TSS_id"))
  
  return(filtered_data)
}





filtered_raw_data <- filter_TSS_sites_by_replicates_new(raw_data, 
                                                       min_count = 50, 
                                                       min_replicates = 5, 
                                                       min_count_all_reps = 10) 
# Retained 26,506 TSS sites passing min_count_all_reps = 10 
# Filtered to 24,072 TSS sites across 15,357 event IDs meeting thresholds.

# --- Prepare Analysis Data Frame ---
TSS_levels_sorted <- unique(raw_data$TSS_id) %>%
  str_sort(numeric = TRUE)

time_levels_sorted <- c('t00', 't04', 't30')
filtered_raw_data <- filtered_raw_data %>%
  mutate(
    timepoint = str_replace(timepoint, "^iPSC", "t00"),
    timepoint = str_replace(timepoint, "^NPC", "t04"),
    timepoint = str_replace(timepoint, "^CN",  "t30")
  )

# Apply the defined factor levels and prepare the main data frame for analysis.
analysis_data <- filtered_raw_data %>%
  mutate(
    BioSampleID = factor(timepoint),
    PolyA_Site = factor(TSS_id, levels = TSS_levels_sorted),
    Time = factor(sub("_.*", "", BioSampleID), levels = time_levels_sorted)
  )

# --- Reshape Data to Wide Format ---
analysis_data_wide <- analysis_data %>%
  pivot_wider(
    names_from = type,      # Create new columns named "Included" and "Skipped"
    values_from = total,    # Populate these columns with values from the 'total' column
    values_fill = 0         # Fill missing combinations with 0
  )

# --- Add Pseudocount and Calculate PSI ---
# A pseudocount of 1 is added to both Included and Skipped reads.
analysis_data_wide <- analysis_data_wide %>%
  mutate(
    Included = Included + 1,
    Skipped = Skipped + 1,
    psi = Included / (Included + Skipped) # Calculate Percent Spliced In (PSI) for each sample
  ) 


write.csv(analysis_data_wide, "./code/AS_APA/output/TSS_counts_wide_new_filter.csv", 
          row.names = F)





# PREPARE EVENTS ----------------------------------------------------------
analysis_data_wide <- read.csv("./code/output/TSS_counts_wide_new_filter.csv")
#Filter to event_ids w/ >1 TSS
filtered_analysis_data_wide <- analysis_data_wide %>%
  group_by(event_id) %>%
  filter(n_distinct(TSS_id) > 1) %>%
  ungroup()
length(unique(filtered_analysis_data_wide$event_id)) #6938
  
#Categorize events to 3 groups (those that pass global only, pass interaction only, and both)

#1) Find events that pass global dPSI:
analysis_data_wide_summary_global <- filtered_analysis_data_wide %>%
  group_by(event_id, TSS_id) %>%
  summarise(mean_psi = mean(psi)) %>% 
  ungroup() %>%
  group_by(event_id) %>%
  mutate(min_psi = min(mean_psi), 
         max_psi = max(mean_psi), 
         dpsi = max_psi - min_psi)
analysis_data_wide_summary_global_filter <- analysis_data_wide_summary_global %>%
  group_by(event_id) %>%
  filter(dpsi >= 0.1)
filtered_analysis_data_wide_global <- filtered_analysis_data_wide %>%
  filter(event_id %in% analysis_data_wide_summary_global_filter$event_id)


#2) Find events that pass interaction dPSI: (test delta for each tp, across TSS site)
analysis_data_wide_summary_int <- filtered_analysis_data_wide %>%
  group_by(event_id, TSS_id, Time) %>%
  summarise(mean_psi = mean(psi)) %>% 
  ungroup() %>%
  group_by(event_id, Time) %>%
  mutate(min_psi = min(mean_psi), 
         max_psi = max(mean_psi), 
         dpsi = max_psi - min_psi)
analysis_data_wide_summary_int_filter <- analysis_data_wide_summary_int %>%
  group_by(event_id, Time) %>%
  filter(dpsi >= 0.1)

#Also make sure there is some variance per polyA site amongst timepoints:
analysis_data_wide_summary_for_int <- filtered_analysis_data_wide %>%
  group_by(event_id, TSS_id, Time) %>%
  summarise(mean_psi = mean(psi)) %>% 
  ungroup() %>%
  group_by(event_id, TSS_id) %>%
  mutate(min_psi = min(mean_psi), 
         max_psi = max(mean_psi), 
         dpsi = max_psi - min_psi)
analysis_data_wide_summary_int_total_filter <- analysis_data_wide_summary_for_int %>%
  group_by(event_id) %>%
  filter(dpsi >= 0.1)

filtered_analysis_data_wide_int_v2 <- filtered_analysis_data_wide %>%
  filter(event_id %in% analysis_data_wide_summary_int_filter$event_id &
           event_id %in% analysis_data_wide_summary_int_total_filter$event_id )


#3) 3 sets of events. 
global_only_events <- setdiff(filtered_analysis_data_wide_global$event_id, filtered_analysis_data_wide_int_v2$event_id) #758
int_only_events <- setdiff(filtered_analysis_data_wide_int_v2$event_id, filtered_analysis_data_wide_global$event_id) #351
both_events <- intersect(filtered_analysis_data_wide_int_v2$event_id, filtered_analysis_data_wide_global$event_id) #1571

all_global <- unique(filtered_analysis_data_wide_global$event_id) #2329
all_int <- unique(filtered_analysis_data_wide_int_v2$event_id) #1922



# QUASI BINOMIAL -----------------------------------------------------
events <- unique(filtered_analysis_data_wide$event_id)
test <- filtered_analysis_data_wide %>% filter(event_id %in% events)

length(unique(test$event_id)) #6938

fit_both_models_and_extract <- function(event_data, event_id) {
  tryCatch({
    # --- Preprocessing ---
    event_data <- event_data %>%
      mutate(
        Time = factor(Time, levels = c("t00", "t04", "t30")),
        TSS = factor(TSS_id, levels = str_sort(unique(TSS_id), numeric = TRUE))
      )
    
    model_family <- "quasibinomial"
    convergence_warning <- FALSE
    
    # --- Fit models and capture convergence warnings ---
    fit_reduced <- withCallingHandlers(
      glm(
        cbind(Included, Skipped) ~ Time,
        data = event_data,
        family = quasibinomial()
      ),
      warning = function(w) {
        if (grepl("algorithm did not converge", conditionMessage(w))) {
          convergence_warning <<- TRUE
        }
      }
    )
    
    fit_full <- withCallingHandlers(
      glm(
        cbind(Included, Skipped) ~ Time + TSS,
        data = event_data,
        family = quasibinomial()
      ),
      warning = function(w) {
        if (grepl("algorithm did not converge", conditionMessage(w))) {
          convergence_warning <<- TRUE
        }
      }
    )
    
    fit_interaction <- withCallingHandlers(
      glm(
        cbind(Included, Skipped) ~ Time * TSS,
        data = event_data,
        family = quasibinomial()
      ),
      warning = function(w) {
        if (grepl("algorithm did not converge", conditionMessage(w))) {
          convergence_warning <<- TRUE
        }
      }
    )
    
    # --- LRTs via F-test ---
    lrt_global <- anova(fit_reduced, fit_full, test = "F")
    lrt_interaction <- anova(fit_full, fit_interaction, test = "F")
    
    global_p <- lrt_global$`Pr(>F)`[2]
    interaction_p <- lrt_interaction$`Pr(>F)`[2]
    
    # --- Convergence flag ---
    model_converged <- isTRUE(fit_interaction$converged) && !convergence_warning
    
    # --- Coefficients from interaction model ---
    coef_df <- summary(fit_interaction)$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("Term") %>%
      mutate(
        event_id = event_id,
        global_LRT_p = global_p,
        interaction_LRT_p = interaction_p,
        interaction_model_converged = model_converged,
        model_family = model_family
      )
    
    # --- Estimated Marginal Means from interaction model ---
    emm_fit <- emmeans(fit_interaction, specs = ~ Time | TSS)
    emm_df <- as.data.frame(emm_fit) %>%
      mutate(
        event_id = event_id,
        global_LRT_p = global_p,
        interaction_LRT_p = interaction_p,
        interaction_model_converged = model_converged,
        model_family = model_family
      )
    
    return(list(coef = coef_df, emmeans = emm_df))
    
  }, error = function(e) {
    warning(paste("Modeling failed for event_id:", event_id, "—", e$message))
    
    coef_fail <- tibble(
      event_id = event_id,
      Term = NA,
      Estimate = NA,
      `Std. Error` = NA,
      `t value` = NA,
      `Pr(>|t|)` = NA,
      global_LRT_p = NA_real_,
      interaction_LRT_p = NA_real_,
      interaction_model_converged = FALSE,
      model_family = "quasibinomial"
    )
    
    emm_fail <- tibble(
      event_id = event_id,
      Time = NA,
      TSS = NA,
      prob = NA,
      SE = NA,
      df = NA,
      lower.CL = NA,
      upper.CL = NA,
      global_LRT_p = NA_real_,
      interaction_LRT_p = NA_real_,
      interaction_model_converged = FALSE,
      model_family = "quasibinomial"
    )
    
    return(list(coef = coef_fail, emmeans = emm_fail))
  })
}



# --- Run on all event_ids ---
split_data <- test %>%
  group_by(event_id) %>%
  group_split()

event_ids <- map_chr(split_data, ~ unique(.x$event_id))

# Run loop and store separately
model_outputs <- map2(split_data, event_ids, fit_both_models_and_extract)

# Extract and bind
results_df <- map_dfr(model_outputs, "coef")
emmeans_df <- map_dfr(model_outputs, "emmeans")

#FDR correction per event (i.e. not per row in the df)
fdr_table_both <- results_df %>%
  group_by(event_id) %>%
  dplyr::slice(1) %>%  # One row per event (we only need one to hold the p-values)
  ungroup() %>%
  mutate(
    global_LRT_FDR = p.adjust(global_LRT_p, method = "fdr"),
    interaction_LRT_FDR = p.adjust(interaction_LRT_p, method = "fdr")
  ) %>%
  dplyr::select(event_id, 
         global_LRT_FDR,
         interaction_LRT_FDR)
# Join back to full coefficient table
all_coefs_with_fdr_both <- results_df %>%
  left_join(fdr_table_both, by = "event_id")


results_pass <- all_coefs_with_fdr_both %>% filter(interaction_model_converged == T)
length(unique(results_pass$event_id)) #6921
write.csv(results_pass, "./code/AS_APA/output/all_events_quasi_TSS.csv", row.names = F)


emmeans_df <- emmeans_df[,-c(9,10)]
emmeans_pass <- emmeans_df %>% filter(interaction_model_converged == T)
length(unique(emmeans_pass$event_id)) #6921
write.csv(emmeans_pass, "./code/AS_APA/output/all_events_emmeans_quasi_TSS.csv", row.names = F)



# ANALYZING EVENTS --------------------------------------------------------
# Apply effect size filter:
both_results <- read.csv("./code/AS_APA/output/all_events_quasi_TSS.csv")

both_results %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #6356
both_results %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #1864
both_results %>% filter(global_LRT_FDR <= 0.05 &
                          interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #1718

both_results_filter <- both_results %>% filter(event_id %in% both_events) #1568
both_results_filter %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #1560
both_results_filter %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #456
both_results_filter %>% filter(global_LRT_FDR <= 0.05 &
                          interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #453


both_results_all_global <- both_results %>% filter(event_id %in% all_global) 
both_results_all_global %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #2304


both_results_all_int <- both_results %>% filter(event_id %in% all_int)
both_results_all_int %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #616


# Compare sig events to sig events from Fisher's - can only compare global
results_list <- read.csv("./code/AS_APA/output/Fisher_result_TSS.csv")

results_list <- results_list %>%
  mutate(
    gene_id = sub(";.*", "", event_id)
    
  )
global_results <- both_results %>% filter(event_id %in% all_global) 
global_results_events <- global_results %>% filter(global_LRT_FDR <= 0.05)
length(unique(global_results_events$event_id)) #2304 sig global events
global_results_events_unique <- unique(global_results_events$event_id)
#   
#Here are the Fisher's events w/ at least one padj < 0.05 and |logOR|>0.05
sig_results #1885 unique
sig_results_unique <- unique(sig_results$event_id)

##But the 'universe' of testable events is different in each:
#Events passing read counts filters for quasi:
quasi_universe <- unique(filtered_analysis_data_wide$event_id) #6938
fisher_universe <- unique(results_list$event_id) #2093

overlapping_universe <- intersect(quasi_universe, fisher_universe) #1885

global_results_events_unique2 <- global_results_events_unique[global_results_events_unique %in% overlapping_universe] #1316
sig_results_unique2 <- sig_results_unique[sig_results_unique %in% overlapping_universe] #1707

length(intersect(global_results_events_unique2, sig_results_unique2)) #1279
length(setdiff(global_results_events_unique2, sig_results_unique2)) #37
length(setdiff(sig_results_unique2, global_results_events_unique2)) #428

#Thus, of the events tested in both analyses (i.e. passing read count filters) = 1885 total
#Then, filtered to sig events
#1279 overlapping (out of 1279+37+428 = 1744 total sig)
#73%



# PLOTTING ----------------------------------------------------
# Load in data
both_results <- read.csv("./code/AS_APA/output/all_events_quasi_TSS.csv")
both_emmeans <- read.csv("./code/AS_APA/output/all_events_emmeans_quasi_TSS.csv")

filtered_df_ES <- read.csv("./code/AS_APA/input/transcript_TSS_sites.csv")
colnames(filtered_df_ES)[5] <- "TSS_id"

counts_annotated <- read.csv("./code/AS_APA/input/tr_TSS_counts.csv")


#All events with at least 1 type of coordination: #unique set of global + int
both_results_all_global <- both_results %>% filter(event_id %in% all_global) 
all_sig_global_events <- both_results_all_global %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)


both_results_all_int <- both_results %>% filter(event_id %in% all_int)
all_sig_int_events <- both_results_all_int %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)

all_sig_events <- unique(union(all_sig_global_events, all_sig_int_events) ) #2467

filtered_analysis_data_wide <- filtered_analysis_data_wide %>% 
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
genes <- unique(filtered_analysis_data_wide$gene_id) #1966

sig_events_df <- filtered_analysis_data_wide %>%
  filter(event_id %in% all_sig_events) %>%
  distinct(event_id, gene_id)
genes_sig <- unique(sig_events_df$gene_id) #1277

#Final event_table:
event_table <- filtered_analysis_data_wide %>% filter(event_id %in% all_sig_events)


# COMPARE TO 3' -----------------------------------------------------------
sig_3 <- read.csv("./code/AS_APA/output/pA_sig_events.csv")

genes_sig_5 <- genes_sig
all_sig_events_5 <- all_sig_events

genes_sig_3 <- unique(sig_3$gene_id)
all_sig_events_3 <- unique(sig_3$event_id)

##after running 3' code:
library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)


gene_list <- list(
  `genes_sig_5` = genes_sig_5,
  `genes_sig_3` = genes_sig_3
)

# genes_venn <- ggVennDiagram(gene_list) +
#   scale_x_continuous(expand = expansion(mult = 0.15))
genes_venn <- venn.diagram(gene_list, category.names = c("", ""),
             filename = NULL,
             scaled = T,
           #  inverted = T,
             rotation.degree = 180,
            # sub = expression(paste(italic("P"), " = 8.56e-288")),
            # sub.pos = c(0.5,0.11),
             fontfamily  = "sans",
             cat.fontfamily = "sans", 
             sub.fontfamily = "sans",
             fill = c("#00A08A", "#9986A5"), 
             alpha = 0.2, 
             col = c("#00A08A", "#9986A5"), 
             cex = 0.8, 
             cat.cex = 1.3, 
          #   cat.pos = c(-20, 20), 
            # cat.dist = c(0.03, 0.07)
            ) 

# test <- ggdraw() +
#   draw_plot(genes_venn,   x = 0.5, y = 0.775, width = 0.5, height = 0.225)  #B
#   
  
#ggsave("./code/AS_APA/genes_venn_diagram.pdf", genes_venn, width = 7, height = 5)

event_list <- list(
  sig_events_5 = all_sig_events_5,
  sig_events_3 = all_sig_events_3
)

events_venn <- venn.diagram(event_list, category.names = c("", ""),
                           filename = NULL,
                           scaled = T,
                           #  inverted = T,
                           rotation.degree = 180,
                           # sub = expression(paste(italic("P"), " = 8.56e-288")),
                           # sub.pos = c(0.5,0.11),
                           fontfamily  = "sans",
                           cat.fontfamily = "sans", 
                           sub.fontfamily = "sans",
                           fill = c("#00A08A", "#9986A5"), 
                           alpha = 0.2, 
                           col = c("#00A08A", "#9986A5"), 
                           cex = 0.8, 
                           cat.cex = 1.3, 
                           #   cat.pos = c(-20, 20), 
                           # cat.dist = c(0.03, 0.07)
) 

venns <- ggdraw() +
  draw_plot(genes_venn,   x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot(events_venn,   x = 0, y = 0, width = 0.5, height = 1) 
venns


svg(paste0("./figures/Figure_5_venns.svg"), height = 4.5, width = 6)
print(venns)
dev.off()




# GENOMIC DISTANCE --------------------------------------------------------
#To get proximal coord per event_id-TSS site:
unique_prox <- filtered_df_ES %>% distinct(gene_id, TSS_id, cluster_proximal_coord)

#To get all event_ids to polyA sites:
event_table <- event_table %>% distinct(event_id, TSS_id) %>%
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
event_table <- merge(event_table, unique_prox, 
                     by.x = c("gene_id", "TSS_id"), by.y = c("gene_id", "TSS_id"))


#Extract exon coords:
coords <- event_table %>%
  # Extract only the coordinate part of event_id: skip everything up to the first colon (after chr...)
  mutate(coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),  # removes up to second colon
         
         # Now extract just the genomic numbers (this avoids chromosome number)
         coord_nums = str_extract_all(coord_part, "\\d+"),
         coord_nums = lapply(coord_nums, as.numeric), 
         strand = str_extract(event_id, ".$") ) %>%
  
  # Calculate most 5' coord based on strand
  rowwise() %>%
  mutate(
    most_5prime_coord = if (strand == "+") min(coord_nums, na.rm = TRUE) else max(coord_nums, na.rm = TRUE),
    distance_to_5prime = abs(cluster_proximal_coord - most_5prime_coord)
  ) %>%
  ungroup()

coords <- coords %>% distinct(event_id, TSS_id, .keep_all = T)


#Filter to the most proximal TSS per event_id:
coords_filtered <- coords %>% 
  mutate(TSS_num = as.integer(sub("TSS", "", TSS_id))   # extract numeric part
  ) %>%
  group_by(event_id) %>%
  slice_max(order_by = TSS_num, with_ties = FALSE) %>%
  ungroup() #2467

median_distance <- median(coords_filtered$distance_to_5prime) #7,531


###ADD ALL OTHER EVENTS:
quasi_universe <- unique(filtered_analysis_data_wide$event_id) #6938
quasi_others <- setdiff(quasi_universe, all_sig_events) #4471

event_table_others <- filtered_analysis_data_wide %>% filter(event_id %in% quasi_others)

event_table_others <- event_table_others %>% distinct(event_id, TSS_id) %>%
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
event_table_others <- merge(event_table_others, unique_prox, 
                            by.x = c("gene_id", "TSS_id"), by.y = c("gene_id", "TSS_id"))


#Extract exon coords:
coords_others <- event_table_others %>%
  # Extract only the coordinate part of event_id: skip everything up to the first colon (after chr...)
  mutate(coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),  # removes up to second colon
         
         # Now extract just the genomic numbers (this avoids chromosome number)
         coord_nums = str_extract_all(coord_part, "\\d+"),
         coord_nums = lapply(coord_nums, as.numeric), 
         strand = str_extract(event_id, ".$") ) %>%
  
  # Calculate most 3' coord based on strand
  rowwise() %>%
  mutate(
    most_5prime_coord = if (strand == "+") min(coord_nums, na.rm = TRUE) else max(coord_nums, na.rm = TRUE),
    distance_to_5prime = abs(cluster_proximal_coord - most_5prime_coord)
  ) %>%
  ungroup()

coords_others <- coords_others %>% distinct(event_id, TSS_id, .keep_all = T)



#Filter to the most proximal pA per event_id:
coords_filtered_others <- coords_others %>% 
  mutate(TSS_num = as.integer(sub("TSS", "", TSS_id))   # extract numeric part
  ) %>%
  group_by(event_id) %>%
  slice_max(order_by = TSS_num, with_ties = FALSE) %>%
  ungroup() #4471

median_distance_others <- median(coords_filtered_others$distance_to_5prime) #10,319

#combine:
coords_filtered$type <- "Sig"
coords_filtered_others$type <- "Others"

coords_all <- bind_rows(coords_filtered, coords_filtered_others)
# 
# 1. Create a summary data frame for the medians
vlines <- coords_all %>%
  group_by(type) %>%
  summarize(med_dist = median(distance_to_5prime, na.rm = TRUE))

ks_res <- ks.test(
  coords_all$distance_to_5prime[coords_all$type == unique(coords_all$type)[1]],
  coords_all$distance_to_5prime[coords_all$type == unique(coords_all$type)[2]]
)

# 2. Format the P-value string for the plot
# This creates a plotmath expression: P = 1.2 x 10^-5
p_label <- paste0("italic(P) == ", gsub("e", " %*%  10^", format(ks_res$p.value, scientific = TRUE, digits = 2)))

# 3. Create the Plot
genomic_distance <- coords_all %>%
  ggplot(aes(x = distance_to_5prime, fill = type, colour = type)) +
  geom_density(alpha = 0.35, linewidth = 0.3, trim = TRUE) +
  scale_fill_manual(values = c("Sig" = "cornflowerblue", "Others" = "firebrick")) +
  scale_colour_manual(values = c("Sig" = "cornflowerblue", "Others" = "firebrick")) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none", # Specify X and Y (0 to 1)
    legend.title = element_blank(),  # Remove legend title
    axis.text = element_text(color = "black"), 
    axis.ticks = element_line(colour = "black"), 
    axis.title.x = element_markdown()
  ) +
  scale_x_log10(
    labels = trans_format("log10", math_format(10^.x)), 
    limits = c(NA, 843520)
  ) +
  # 1) Thicker vlines
  geom_vline(
    data = vlines, 
    aes(xintercept = med_dist, colour = type), 
    linetype = "dashed",
    linewidth = 0.3, # Increased thickness
    show.legend = FALSE
  ) +
  # 3) Formatted P-value
  annotate(
    "text", 
    x = 12, y = Inf, 
    label = p_label, 
    parse = TRUE,    # This tells ggplot to render the math expression
    hjust = 0.5, vjust = 2, 
    size = 3
  ) +
  labs(
    x = "<span style='font-size:9pt'>Genomic Distance (nt)</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
    y = "Density"
  )

genomic_distance


# genomic_distance <- coords_all %>%
#   ggplot(aes(x = distance_to_3prime, group = type )) +
#   geom_density(alpha = 0.35, linewidth = 0.6, trim = T, fill = "cornflowerblue", colour = "darkblue") +
#   theme_classic(base_size = 9) +
#   theme(legend.position = "none", 
#         axis.text = element_text(color = "black"), 
#         axis.ticks = element_line(colour = "black"), 
#         axis.title.x = element_markdown() ) +
#   scale_x_log10(
#     labels = trans_format("log10", math_format(10^.x)), 
#     limits = c(9, NA)
#   ) +
#   geom_vline(xintercept = median_distance, colour = "black", linetype = "dashed") +
#   labs(
#     x = "<span style='font-size:9pt'>Genomic Distance (nt)</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
#     y = "Density")
# 
# genomic_distance
#  



# NUMBER OF EXONS ---------------------------------------------------------
ORFanage_replaced <- rtracklayer::import("./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_gtf.gtf") #only using 158,844 here

#Filter to only transcripts in the significant event-pA IDs (only for the most proximal pA)
counts_filter <- counts_annotated %>% 
  filter(paste(event_id, TSS_id) %in% paste(coords_filtered$event_id, coords_filtered$TSS_id)) #2467
counts_filter <- counts_filter %>%
  distinct(transcript_id, event_id, .keep_all = T)

#Sum across all samples:
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

tr_count <- tr_count %>%
  mutate(sum_count = rowSums(dplyr::select(., iPSC_1:CN_3_2)))

counts_filter <- merge(counts_filter, tr_count[,c("isoform", "sum_count")], 
                       by.x = "transcript_id", by.y = "isoform")
counts_max <- counts_filter %>%
  group_by(event_id) %>%
  slice_max(sum_count, n = 1, with_ties = FALSE) %>%
  ungroup() #2467


#Use gtf to get number of exons from ... to transcript_end
counts_max <- counts_max %>%
  # Extract only the coordinate part of event_id: skip everything up to the first colon (after chr...)
  mutate(coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),  # removes up to second colon
         # Now extract just the genomic numbers (this avoids chromosome number)
         coord_nums = str_extract_all(coord_part, "\\d+"),
         coord_nums = lapply(coord_nums, as.numeric),
         strand = str_extract(event_id, "[+-]$"),
         chr = str_extract(event_id, "chr([0-9]{1,2}|X|Y)"),
         # Pick number depending on strand
         selected_coord = map2_dbl(coord_nums, strand, ~ if(.y == "+") .x[2] else .x[3]) )


#Now get # exons b/w:
exons <- ORFanage_replaced[ORFanage_replaced$type == "exon"]

# Step 1: create regions GRanges for each event_id
# Make sure start <= end even for negative strand
regions <- GRanges(
  seqnames = as.character(counts_max$chr),
  ranges   = IRanges(
    start = pmin(counts_max$selected_coord, counts_max$transcript_start),
    end   = pmax(counts_max$selected_coord, counts_max$transcript_start)
  ),
  strand = as.character(counts_max$strand),
  transcript_id = counts_max$transcript_id,
  event_id      = counts_max$event_id
)

# Step 2: count exons per region (event_id)
#Include exons that overlap transcript_start (left boundary).
#Exclude exons that overlap selected_coord (right boundary).

num_exons_in_region <- sapply(seq_along(regions), function(i) {
  region_i <- regions[i]
  tx_id    <- mcols(region_i)$transcript_id
  
  # subset exons to this transcript
  exons_tx <- exons[mcols(exons)$transcript_id == tx_id]
  
  # adjust region depending on strand
  if (as.character(strand(region_i)) == "+") {
    region_i_adj <- GRanges(
      seqnames = seqnames(region_i),
      ranges   = IRanges(
        start = start(region_i), 
        end   = end(region_i) -1,
      ),
      strand   = strand(region_i)
    )
  } else {
    region_i_adj <- GRanges(
      seqnames = seqnames(region_i),
      ranges   = IRanges(
        start = start(region_i) + 1, 
        end   = end(region_i)
      ),
      strand   = strand(region_i)
    )
  }
  
  # count overlaps
  length(findOverlaps(region_i_adj, exons_tx, ignore.strand = FALSE))
})

# Step 3: assign back to counts_max
counts_max$num_exons_in_region <- num_exons_in_region

median_num_exon <- median(counts_max$num_exons_in_region) #3


####DO FOR UNCOORD:
counts_filter_others <- counts_annotated %>% 
  filter(paste(event_id, TSS_id) %in% paste(coords_filtered_others$event_id, coords_filtered_others$TSS_id)) 
counts_filter_others <- counts_filter_others %>%
  distinct(transcript_id, event_id, .keep_all = T)

#Sum across all samples:
counts_filter_others <- merge(counts_filter_others, tr_count[,c("isoform", "sum_count")], 
                              by.x = "transcript_id", by.y = "isoform")
counts_max_others <- counts_filter_others %>%
  group_by(event_id) %>%
  slice_max(sum_count, n = 1, with_ties = FALSE) %>%
  ungroup() #4471


#Use gtf to get number of exons from ... to transcript_start
counts_max_others <- counts_max_others %>%
  # Extract only the coordinate part of event_id: skip everything up to the first colon (after chr...)
  mutate(coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),  # removes up to second colon
         # Now extract just the genomic numbers (this avoids chromosome number)
         coord_nums = str_extract_all(coord_part, "\\d+"),
         coord_nums = lapply(coord_nums, as.numeric),
         strand = str_extract(event_id, "[+-]$"),
         chr = str_extract(event_id, "chr([0-9]{1,2}|X|Y)"),
         # Pick number depending on strand
         selected_coord = map2_dbl(coord_nums, strand, ~ if(.y == "+") .x[2] else .x[3]) )



# Step 1: create regions GRanges for each event_id
# Make sure start <= end even for negative strand
regions_others <- GRanges(
  seqnames = as.character(counts_max_others$chr),
  ranges   = IRanges(
    start = pmin(counts_max_others$selected_coord, counts_max_others$transcript_start),
    end   = pmax(counts_max_others$selected_coord, counts_max_others$transcript_start)
  ),
  strand = as.character(counts_max_others$strand),
  transcript_id = counts_max_others$transcript_id,
  event_id      = counts_max_others$event_id
)

# Step 2: count exons per region (event_id)
#Include exons that overlap transcript_start (right boundary).
#Exclude exons that overlap selected_coord (left boundary).

num_exons_in_region_others <- sapply(seq_along(regions_others), function(i) {
  region_i <- regions_others[i]
  tx_id    <- mcols(region_i)$transcript_id
  
  # subset exons to this transcript
  exons_tx <- exons[mcols(exons)$transcript_id == tx_id]
  
  # adjust region depending on strand
  if (as.character(strand(region_i)) == "+") {
    region_i_adj <- GRanges(
      seqnames = seqnames(region_i),
      ranges   = IRanges(
        start = start(region_i), 
        end   = end(region_i) -1,
      ),
      strand   = strand(region_i)
    )
  } else {
    region_i_adj <- GRanges(
      seqnames = seqnames(region_i),
      ranges   = IRanges(
        start = start(region_i) + 1, 
        end   = end(region_i)
      ),
      strand   = strand(region_i)
    )
  }
  
  # count overlaps
  length(findOverlaps(region_i_adj, exons_tx, ignore.strand = FALSE))
})

# Step 3: assign back to counts_max
counts_max_others$num_exons_in_region <- num_exons_in_region_others

median_num_exon_others <- median(counts_max$num_exons_in_region) #3



# num_exons <- counts_max %>%  
#   ggplot(aes(x = num_exons_in_region)) +
#   geom_bar(width = 0.7, alpha = 0.35, linewidth = 0.25,
#            fill = "cornflowerblue", color = "darkblue") +
#   theme_classic(base_size = 9) +
#   coord_cartesian(ylim = c(0, 375)) +
#   labs(
#     x = "<span style='font-size:9pt'>Number of Exons</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
#     y = "# of AS events") +
#   theme(axis.text.x = element_text(color = "black"), 
#         axis.title.x = element_markdown(),
#         axis.text.y = element_text(colour = "black"), 
#         axis.ticks = element_line(colour = "black") ) +
#   geom_vline(xintercept = median_num_exon, colour = "black", linetype = "dashed") 
# 
# num_exons

#combine:
counts_max$type <- "Sig"
counts_max_others$type <- "Others"
counts_all <- bind_rows(counts_max, counts_max_others)

vlines_exons <- counts_all %>%
  group_by(type) %>%
  summarize(med_dist = median(median_num_exon, na.rm = TRUE))
vlines_offset <- vlines_exons %>%
  mutate(med_dist = ifelse(row_number() == 1, 
                           med_dist * 0.96,  # Nudge first line slightly left
                           med_dist * 1.04)) # Nudge second line slightly right

num_exons <- counts_all %>%  
  ggplot(aes(x = num_exons_in_region, fill = type, colour = type)) +
  geom_bar(width = 0.7, alpha = 0.35, linewidth = 0.25) +
  scale_fill_manual(values = c("Sig" = "cornflowerblue", "Others" = "firebrick")) +
  scale_colour_manual(values = c("Sig" = "cornflowerblue", "Others" = "firebrick")) +
  theme_classic(base_size = 9) +
  # coord_cartesian(ylim = c(0, 375)) +
  labs(
    x = "<span style='font-size:9pt'>Number of Exons</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
    y = "# of AS events") +
  theme( legend.key.size = unit(0.30, 'cm'),
         legend.position = c(0.9, 0.8), # Specify X and Y (0 to 1)
         legend.title = element_blank(),  # Remove legend title
         axis.text.x = element_text(color = "black"), 
         axis.title.x = element_markdown(),
         axis.text.y = element_text(colour = "black"), 
         axis.ticks = element_line(colour = "black") ) +
  geom_vline(
    data = vlines_offset, 
    aes(xintercept = med_dist, colour = type), 
    linetype = "dashed",
    linewidth = 0.3, # Increased thickness
    show.legend = FALSE
  ) 

num_exons



# COORDINATED EXON TR REGION ----------------------------------------------
###Are the coordinated exons in 5'UTR, ORF or 3' UTR?
##1) Take the most abundant transcript (WITH THE EXON INCLUDED) from each event (from any pA)

max_tr <- merge(counts_annotated, tr_count[,c("isoform", "sum_count")], 
                by.x = "transcript_id", by.y = "isoform")
max_tr <- max_tr %>% filter(event_id %in% all_sig_events)
max_tr <- max_tr %>%
  group_by(event_id) %>%
  filter(type == "Included") %>%
  slice_max(sum_count, n = 1, with_ties = FALSE) %>%
  ungroup() #2647


#First, extract the exon coordinates:
max_tr <- max_tr %>% 
  mutate(
    # drop everything up to the second colon (keeps coord1-coord2:coord3-coord4:-)
    coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),
    
    # extract just the numbers
    coord_nums = str_extract_all(coord_part, "\\d+") |> lapply(as.numeric),
    
    # strand and chr
    strand = str_extract(event_id, "[+-]$"),
    chr    = str_extract(event_id, "chr([0-9]{1,2}|X|Y|M)")
  ) %>%
  # unpack into separate start/end columns
  mutate(
    exon_start   = sapply(coord_nums, `[`, 2),
    exon_end = sapply(coord_nums, `[`, 3),
  )



# exons_gr: GRanges of exons with transcript_id
exons_gr <- GRanges(
  seqnames = as.character(max_tr$chr),
  ranges   = IRanges(
    start = pmin(max_tr$exon_start, max_tr$exon_end),
    end   = pmax(max_tr$exon_start, max_tr$exon_end)
  ),
  strand = as.character(max_tr$strand),
  transcript_id = max_tr$transcript_id,
  event_id      = max_tr$event_id
)

# features_gr: GRanges of CDS, UTR, etc.
gtf_sub <- ORFanage_replaced[ORFanage_replaced$transcript_id %in% unique(max_tr$transcript_id)] #should be only 2073 transcripts

five_prime_utr_gr <- gtf_sub[gtf_sub$type == "five_prime_utr"]
CDS_gr           <- gtf_sub[gtf_sub$type == "CDS"]
three_prime_utr_gr <- gtf_sub[gtf_sub$type == "three_prime_utr"]

#Split by tr ID:
five_prime_utr_by_tx <- split(five_prime_utr_gr, five_prime_utr_gr$transcript_id)
CDS_by_tx            <- split(CDS_gr, CDS_gr$transcript_id)
three_prime_utr_by_tx <- split(three_prime_utr_gr, three_prime_utr_gr$transcript_id)

label_exon <- function(exon, tx_id) {
  feats <- c()
  
  if (tx_id %in% names(five_prime_utr_by_tx)) {
    f <- five_prime_utr_by_tx[[tx_id]]
    ov <- findOverlaps(exon, f)
    if (length(ov) > 0) feats <- c(feats, "five_prime_UTR")
  }
  
  if (tx_id %in% names(CDS_by_tx)) {
    f <- CDS_by_tx[[tx_id]]
    ov <- findOverlaps(exon, f)
    if (length(ov) > 0) feats <- c(feats, "CDS")
  }
  
  if (tx_id %in% names(three_prime_utr_by_tx)) {
    f <- three_prime_utr_by_tx[[tx_id]]
    ov <- findOverlaps(exon, f)
    if (length(ov) > 0) feats <- c(feats, "three_prime_UTR")
  }
  
  if (length(feats) == 0) return("unlabeled")
  if (length(feats) == 1) return(paste0("within_", feats))
  return(paste("overlaps", paste(feats, collapse="+")))
}

exons_gr$label <- sapply(seq_along(exons_gr), function(i) {
  label_exon(exons_gr[i], exons_gr$transcript_id[i])
})

#
label_counts <- as.data.frame(exons_gr) %>%
  mutate(label = case_when(
    label == "overlaps CDS+three_prime_UTR" ~ "ORF + 3'UTR",
    label == "overlaps five_prime_UTR+CDS" ~ "5'UTR + ORF",
    label == "within_CDS" ~ "ORF",
    label == "within_five_prime_UTR" ~ "5'UTR",
    label == "within_three_prime_UTR" ~ "3'UTR"
  ))
preferred_order <- c("5'UTR", "5'UTR + ORF", "ORF", "ORF + 3'UTR", "3'UTR")
label_counts$label <- factor(label_counts$label, levels = preferred_order)


# Plot
overlap_type_bar <- ggplot(label_counts, aes(x = label, fill = label)) +
  geom_bar() +  
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(hjust = 0.5, color = "black"), 
        axis.text.y = element_text(colour = "black"), 
        axis.ticks = element_line(colour = "black"), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 9, face = "bold")) +
  scale_x_discrete(labels = label_wrap(width = 6)) +
  labs(y = "# of AS events", x = "Transcript Region",
       title = "Location of Coordinated Exons")
overlap_type_bar



# FIND GOOD EXAMPLE -------------------------------------------------------
SFARI <- read.csv("./data/SFARI-Gene_genes_04-03-2025release_04-15-2025export.csv")
SFARI_1 <- SFARI %>% filter(gene.score == 1)
  
event_table_SFARI <- event_table %>% filter(gene_id %in% SFARI$gene.symbol )
event_table_SFARI_1 <- event_table %>% filter(gene_id %in% SFARI_1$gene.symbol )

plotting_events <- unique(event_table_SFARI_1$event_id)
plotting_events <- unique(event_table_SFARI$event_id)

pa_levels_sorted <- unique(raw_data$TSS_id) %>%
  str_sort(numeric = TRUE)

time_levels_sorted <- c('t00', 't04', 't30')

##PSI plots for visualization:
plot_data <- filtered_analysis_data_wide %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("t00", "t04", "t30"))),
    PolyA_Site = factor(as.factor(PolyA_Site), levels = rev(pa_levels_sorted))
  )
# Calculate mean PSI for each group to determine the height of the bars.
summary_data <- filtered_analysis_data_wide %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(mean_psi = mean(psi, na.rm = TRUE), .groups = 'drop') %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("t00", "t04", "t30"))),
    PolyA_Site = factor(as.factor(PolyA_Site), levels = rev(pa_levels_sorted))
  )

both_results <- both_results %>% 
  mutate(gene_id = str_extract(event_id, "^[^;]+"))

facet_labels <- both_results %>%
  group_by(event_id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      gene_id, "\n",
      "Global FDR = ", signif(global_LRT_FDR, 3), "\n",
      "Interaction FDR = ", signif(interaction_LRT_FDR, 3)
    )
  ) %>%
  dplyr::select(event_id, label)
label_vec <- setNames(facet_labels$label, facet_labels$event_id)


##Add SD:
summary_stats <- plot_data %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(
    mean_psi = mean(psi, na.rm = TRUE),
    sd_psi = sd(psi, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()


custom_colors <- c(
  "t00" = "#D4C1EC",
  "t04"  = "#8E8EEA",
  "t30"   = "#5F58EA")

# Generate the PSI plot.
psi_plot <- summary_stats %>% filter(event_id %in% plotting_events) %>%
  ggplot(aes(y = PolyA_Site, fill = Time)) +
  # Bars represent the mean PSI value for each group.
  geom_col(
    aes(x = mean_psi),
    position = position_dodge(width = 0.9),
    alpha = 0.7 # Use transparency to make overlaid points more visible.
  ) +
  geom_errorbar(aes(xmin = mean_psi - sd_psi, xmax = mean_psi + sd_psi), 
                position = position_dodge(width = 0.9),  # <- this fixes the alignment
                width = 0.2) +
  
  # Points represent the PSI value for each individual biological replicate.
  geom_jitter(
    data = plot_data %>% filter(event_id %in% plotting_events),
    aes(x = psi),
    position = position_jitterdodge(
      jitter.width = 0.2, # Control horizontal spread of points.
      dodge.width = 0.9   # Ensure points align with their corresponding bars.
    ),
    shape = 21, color = "black", size = 1
  ) +
  coord_cartesian(xlim = c(0, 1)) +  # <- keeps all error bars visible
  scale_fill_manual(values = custom_colors, 
                    breaks = c("t00", "t04", "t30")) +
  labs(x = "PSI",
       y = ""
  ) +
  theme_classic(base_size = 9) +
  facet_wrap(~event_id, scales = "free", labeller = as_labeller(label_vec), 
             ncol = 3) +
  theme( legend.key.size = unit(0.30, 'cm'),
         legend.position = c(0.85, 0.32),
         legend.title = element_blank()     ) 
psi_plot



pdf(paste0("./figures/tss_sfari_1.pdf"), height = 45, width = 8)
print(psi_plot)
dev.off()


pdf(paste0("./figures/tss_sfari.pdf"), height = 100, width = 8)
print(psi_plot)
dev.off()

##good global example: TSC2;SE:chr16:2050486-2053342:2053452-2054296:+
##good interaction example: SYNCRIP;SE:chr6:85640564-85641292:85641451-85642797:-


# PLOTTING TSC2 -----------------------------------------------------------
TSC2 <- counts_annotated %>% filter(event_id == "TSC2;SE:chr16:2050486-2053342:2053452-2054296:+") %>%
  distinct(event_id, gene_id, transcript_id, TSS_id) %>%
  filter(TSS_id %in% c("TSS1", "TSS2"))

counts_test <- tr_count %>% filter(isoform %in% TSC2$transcript_id)
  

#Plot 3 sep plots, 1 for each timepoint:
timepoints <- c("sum_t00", "sum_t04", "sum_t30")
TSS_ids <- TSC2 %>% distinct(TSS_id) %>% pull(TSS_id)
desired_order <- c("TSS1", "TSS2")
TSS_ids <- desired_order[desired_order %in% TSS_ids]


TSC2_gtf_1 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.91603.1") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(TSC2_gtf_1, "./code/AS_APA/input/TSC2_gtf_1.gtf")

TSC2_gtf_2 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.91603.2") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(TSC2_gtf_2, "./code/AS_APA/input/TSC2_gtf_2.gtf")


txdb_1 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/TSC2_gtf_1.gtf", format = "gtf")
txdb_2 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/TSC2_gtf_2.gtf", format = "gtf")


fill_colors <- rep("steelblue", 11) 
fill_colors[c(5)] <- "red" ##This is manual...

gviz_fontsize <- 7
# 

txTr_1 <- GeneRegionTrack(txdb_1,
                            col = NA,
                            fill = fill_colors,
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "TSS1")
displayPars(txTr_1) <- list(transcript = "TSS1",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_2 <- GeneRegionTrack(txdb_2,
                            col = NA,
                            fill = "steelblue",
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "TSS2")
displayPars(txTr_2) <- list(transcript = "TSS2",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 


# Combine tracks: Genome axis, data tracks, and transcripts

# Plot
plotTracks(c(GenomeAxisTrack(name = ""), txTr_1, txTr_2),
           main = "TSC2",
           from = 2047356, 
           to = 2059072,
           showTitle = T,
           showId = T,
           cex.main = 0.8,
           just.group = "left",
           background.title = "white",
           col.title        = "black",
           background.panel = "white",
           col.axis         = "black",
           col.frame        = "white" )

# Save Gviz as grob
gviz_grob <- grid::grid.grabExpr(
  plotTracks(c(GenomeAxisTrack(name = "", scale = 0.25), #DataTracks,
               txTr_1, txTr_2),
             # main = "TSC2 at t00",
             from = 2047356, 
             to = 2059072,
             showTitle = T,
             cex.main = 0.8,
             just.group = "left",
             background.title = "white",
             col.title        = "black",
             background.panel = "white",
             col.axis         = "white", #hack to hide the axis
             col.frame        = "white", 
             sizes = c(0.1, 0.1, 0.1), #0.925
             title.width = 0.5)
)



# PLOTTING SYNCRIP ---------------------------------------------------------
#Representative tr is the most highly expressed:
SYNCRIP <-  counts_annotated %>% filter(event_id == "SYNCRIP;SE:chr6:85640564-85641292:85641451-85642797:-") %>%
  distinct(event_id, gene_id, transcript_id, TSS_id) %>%
  filter(TSS_id %in% c("TSS1", "TSS2", "TSS3", "TSS4"))

counts_test2 <- tr_count %>% filter(isoform %in% SYNCRIP$transcript_id)
SYNCRIP <- merge(SYNCRIP, counts_test2[,c('isoform', 'sum_count')], 
                 by.x = "transcript_id", by.y = "isoform")
#Representative tr is the most highly expressed:
# TSS1: PB.44678.1077
#TSS2: PB.44678.1124
#TSS3: PB.44678.1239
#TSS4: PB.44678.1005


SYNCRIP_gtf_1077 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.44678.1077") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SYNCRIP_gtf_1077, "./code/AS_APA/input/SYNCRIP_gtf_1077.gtf")

SYNCRIP_gtf_1124 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.44678.1124") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SYNCRIP_gtf_1124, "./code/AS_APA/input/SYNCRIP_gtf_1124.gtf")


SYNCRIP_gtf_1239 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.44678.1239") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SYNCRIP_gtf_1239, "./code/AS_APA/input/SYNCRIP_gtf_1239.gtf")


SYNCRIP_gtf_1005 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.44678.1005") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SYNCRIP_gtf_1005, "./code/AS_APA/input/SYNCRIP_gtf_1005.gtf")


txdb_1077 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/SYNCRIP_gtf_1077.gtf", format = "gtf")
txdb_1124 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/SYNCRIP_gtf_1124.gtf", format = "gtf")
txdb_1239 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/SYNCRIP_gtf_1239.gtf", format = "gtf")
txdb_1005 <- txdbmaker::makeTxDbFromGFF("./code/AS_APA/input/SYNCRIP_gtf_1005.gtf", format = "gtf")



fill_colors <- rep("steelblue", 15) 
fill_colors[c(12)] <- "red" ##This is manual...

gviz_fontsize <- 7

txTr_1077 <- GeneRegionTrack(txdb_1077,
                            col = NA,
                            fill = fill_colors,
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "TSS1")
displayPars(txTr_1077) <- list(transcript = "TSS1",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_1124 <- GeneRegionTrack(txdb_1124,
                            col = NA,
                            fill = "steelblue",
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "TSS2"
)
displayPars(txTr_1124) <- list(transcript = "TSS2",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_1239 <- GeneRegionTrack(txdb_1239,
                             col = NA,
                             fill = fill_colors,
                             fontcolor.group = "transparent",
                             showId = F,
                             collapse = FALSE,
                             min.distance = 0,
                             name = "TSS3"
)
displayPars(txTr_1239) <- list(transcript = "TSS3",
                               showId = FALSE, 
                               fontsize = gviz_fontsize, 
                               min.height = 25*0.3) 

txTr_1005 <- GeneRegionTrack(txdb_1005,
                             col = NA,
                             fill = fill_colors,
                             fontcolor.group = "transparent",
                             showId = F,
                             collapse = FALSE,
                             min.distance = 0,
                             name = "TSS4"
)
displayPars(txTr_1005) <- list(transcript = "TSS4",
                               showId = FALSE, 
                               fontsize = gviz_fontsize, 
                               min.height = 25*0.3) 




# Plot
plotTracks(c(GenomeAxisTrack(name = ""),
             txTr_1077, txTr_1124, txTr_1239, txTr_1005),
           from = 85636834, 
           to = 85643981,
           reverseStrand = T,
           main = "SYNCRIP",
           showTitle = T,
           showId = T,
           cex.main = 0.8,
           just.group = "left",
           background.title = "white",
           col.title        = "black",
           background.panel = "white",
           col.axis         = "black",
           col.frame        = "white" )

# Save Gviz as grob
gviz_grob_SYNCRIP <- grid::grid.grabExpr(
  plotTracks(c(GenomeAxisTrack(name = "", scale = 0.25), 
               txTr_1077, txTr_1124, txTr_1239, txTr_1005),
             # main = "SYNCRIP at t00",
             from = 85636834, 
             to = 85643981,
             reverseStrand = T,
             showTitle = T,
             cex.main = 0.8,
             just.group = "left",
             background.title = "white",
             col.title        = "black",
             background.panel = "white",
             col.axis         = "white", #hack to hide the axis
             col.frame        = "white", 
             sizes = c(0.1, 0.1, 0.1, 0.1, 0.1), #0.925
             title.width = 0.5)
)


# PLOTTING TSC2 & SYNCRIP PSI -----------------------------------------------
plotting_events <- c("TSC2;SE:chr16:2050486-2053342:2053452-2054296:+", 
                     "SYNCRIP;SE:chr6:85640564-85641292:85641451-85642797:-")
pa_levels_sorted <- unique(raw_data$TSS_id) %>%
  str_sort(numeric = TRUE)

time_levels_sorted <- c('t00', 't04', 't30')

##PSI plots for visualization:
plot_data <- filtered_analysis_data_wide %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("t00", "t04", "t30"))),
    PolyA_Site = factor(as.factor(PolyA_Site), levels = rev(pa_levels_sorted))
  )
# Calculate mean PSI for each group to determine the height of the bars.
summary_data <- filtered_analysis_data_wide %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(mean_psi = mean(psi, na.rm = TRUE), .groups = 'drop') %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("t00", "t04", "t30"))),
    PolyA_Site = factor(as.factor(PolyA_Site), levels = rev(pa_levels_sorted))
  )

both_results <- both_results %>% 
  mutate(gene_id = str_extract(event_id, "^[^;]+"))

facet_labels <- both_results %>%
  group_by(event_id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      gene_id, "\n",
      "Global FDR = ", signif(global_LRT_FDR, 3), "\n",
      "Interaction FDR = ", signif(interaction_LRT_FDR, 3)
    )
  ) %>%
  dplyr::select(event_id, label)
label_vec <- setNames(facet_labels$label, facet_labels$event_id)


##Add SD:
summary_stats <- plot_data %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(
    mean_psi = mean(psi, na.rm = TRUE),
    sd_psi = sd(psi, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()


custom_colors <- c(
  "t00" = "#D4C1EC",
  "t04"  = "#8E8EEA",
  "t30"   = "#5F58EA")

# Generate the PSI plot.
psi_plot <- summary_stats %>% filter(event_id %in% plotting_events) %>%
  ggplot(aes(y = PolyA_Site, fill = Time)) +
  # Bars represent the mean PSI value for each group.
  geom_col(
    aes(x = mean_psi),
    position = position_dodge(width = 0.9),
    alpha = 0.7 # Use transparency to make overlaid points more visible.
  ) +
  geom_errorbar(aes(xmin = mean_psi - sd_psi, xmax = mean_psi + sd_psi), 
                position = position_dodge(width = 0.9),  # <- this fixes the alignment
                width = 0.2) +
  
  # Points represent the PSI value for each individual biological replicate.
  geom_jitter(
    data = plot_data %>% filter(event_id %in% plotting_events),
    aes(x = psi),
    position = position_jitterdodge(
      jitter.width = 0.2, # Control horizontal spread of points.
      dodge.width = 0.9   # Ensure points align with their corresponding bars.
    ),
    shape = 21, color = "black", size = 1
  ) +
  coord_cartesian(xlim = c(0, 1)) +  # <- keeps all error bars visible
  scale_fill_manual(values = custom_colors, 
                    breaks = c("t00", "t04", "t30")) +
  labs(x = "PSI",
       y = ""
  ) +
  theme_classic(base_size = 9) +
  facet_wrap(~event_id, scales = "free", labeller = as_labeller(label_vec), 
             ncol = 1) +
  theme( legend.key.size = unit(0.30, 'cm'),
         legend.position = c(0.85, 0.32),
         legend.title = element_blank()     ) 
psi_plot



# PLOTTING AGO1 + SMARCA4 COUNTS ------------------------------------------
events_pas <- paste0(event_table$event_id, event_table$TSS_id) %>% unique()


raw_data2 <- raw_data %>%
  mutate(id = paste0(event_id, TSS_id)) %>%
  filter(id %in% events_pas)  %>%
  mutate(
    timepoint = str_replace(timepoint, "^iPSC", "t00"),
    timepoint = str_replace(timepoint, "^NPC", "t04"),
    timepoint = str_replace(timepoint, "^CN",  "t30")
  )

raw_data2 <- raw_data2 %>%
  group_by(event_id, TSS_id, timepoint) %>%
  mutate(total_counts = sum(total))
raw_data2 <- raw_data2 %>% distinct(event_id, TSS_id, timepoint, .keep_all = T)
raw_data2$PolyA_Site <- raw_data2$TSS_id

counts_data <- raw_data2 %>%
  mutate(
    BioSampleID = factor(timepoint),
    PolyA_Site = factor(TSS_id, levels = rev(pa_levels_sorted)),
    Time = factor(sub("_.*", "", BioSampleID), levels = rev(time_levels_sorted))
  )
# Calculate mean counts for each group to determine the height of the bars.
#counts_data_filtered <- counts_data %>% filter(event_id %in% plotting_events)
summary_data_counts <- counts_data %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(mean_counts = mean(total_counts, na.rm = TRUE), .groups = 'drop') %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("t00", "t04", "t30"))),
    PolyA_Site = factor(as.factor(PolyA_Site), levels = rev(pa_levels_sorted))
  )


##Add SD:
summary_stats_counts <- counts_data %>%
  group_by(event_id, PolyA_Site, Time) %>%
  summarise(
    mean_counts = mean(total_counts, na.rm = TRUE),
    sd_counts = sd(total_counts, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()


custom_colors <- c(
  "t00" = "#D4C1EC",
  "t04"  = "#8E8EEA",
  "t30"   = "#5F58EA")

# Generate the counts plot.
counts_plot <- summary_stats_counts %>% filter(event_id %in% plotting_events) %>%
  ggplot(aes(y = PolyA_Site, fill = Time)) +
  # Bars represent the mean counts value for each group.
  geom_col(
    aes(x = mean_counts),
    position = position_dodge(width = 0.9),
    alpha = 0.7 # Use transparency to make overlaid points more visible.
  ) +
  geom_errorbar(aes(xmin = mean_counts - sd_counts, xmax = mean_counts + sd_counts), 
                position = position_dodge(width = 0.9),  # <- this fixes the alignment
                width = 0.2) +
  
  # Points represent the counts value for each individual biological replicate.
  geom_jitter(
    data = counts_data %>% filter(event_id %in% plotting_events),
    aes(x = total_counts),
    position = position_jitterdodge(
      jitter.width = 0.2, # Control horizontal spread of points.
      dodge.width = 0.9   # Ensure points align with their corresponding bars.
    ),
    shape = 21, color = "black", size = 1
  ) +
  # coord_cartesian(xlim = c(0, 1)) +  # <- keeps all error bars visible
  scale_fill_manual(values = custom_colors, 
                    breaks = c("t00", "t04", "t30")) +
  labs(x = "counts",
       y = ""
  ) +
  theme_classic(base_size = 9) +
  facet_wrap(~event_id, scales = "free", labeller = as_labeller(label_vec), 
             ncol = 1) +
  theme( legend.key.size = unit(0.30, 'cm'),
         legend.position = c(0.85, 0.32),
         legend.title = element_blank()     ) 
counts_plot


###COMBINED PLOT:
combined_summary <- bind_rows(
  summary_stats %>%
    mutate(
      value = mean_psi,
      sd = sd_psi,
      metric = "PSI"
    ) %>%
    dplyr::select(event_id, PolyA_Site, Time, value, sd, metric),
  
  summary_stats_counts %>%
    mutate(
      value = mean_counts,
      sd = sd_counts,
      metric = "Counts"
    ) %>%
    dplyr::select(event_id, PolyA_Site, Time, value, sd, metric)
)

combined_points <- bind_rows(
  plot_data %>%
    dplyr::select(event_id, PolyA_Site, Time, psi) %>%
    dplyr::rename(value = psi) %>%
    dplyr::mutate(metric = "PSI"),
  
  counts_data %>%
    dplyr::select(event_id, PolyA_Site, Time, total_counts) %>%
    dplyr::rename(value = total_counts) %>%
    dplyr::mutate(metric = "Counts")
)

combined_summary <- combined_summary %>%
  mutate(gene_id = stringr::str_extract(event_id, "^[^;]+"))

combined_points <- combined_points %>%
  mutate(gene_id = stringr::str_extract(event_id, "^[^;]+"))

# 
# format_sci_expr <- function(x, digits = 2) {
#   if (is.na(x)) return(NA)
#   if (x == 0) return(0)
#   
#   exp <- floor(log10(x))
#   base <- signif(x / 10^exp, digits)
#   
#   bquote(.(base) %*% 10^.(exp))
# }
# 
# facet_labels <- both_results %>%
#   group_by(event_id) %>%
#   dplyr::slice(1) %>%
#   ungroup() %>%
#   mutate(
#     label = paste0(
#       gene_id, "\n",
#       "Global FDR = ", format_sci(global_LRT_FDR),
#       "Interaction FDR = ", format_sci(interaction_LRT_FDR)
#       
#     )
#   ) %>%
#   dplyr::select(event_id, label)
# label_vec <- setNames(facet_labels$label, facet_labels$event_id)
# 
# 
# 
# facet_labels2 <- facet_labels %>%
#   tidyr::crossing(metric = c("PSI", "Counts")) %>%
#   mutate(label = paste0(label, "\n", metric))
# 
# label_vec2 <- setNames(facet_labels2$label,
#                        paste(facet_labels2$event_id, facet_labels2$metric))
# 


plot_event <- function(event) {
  
  # ---------- subset data ----------
  sum_sub <- combined_summary %>% filter(event_id == event)
  pts_sub <- combined_points %>% filter(event_id == event)
  
  # ---------- PSI plot (keeps y-axis labels) ----------
  psi_plot <- ggplot(
    sum_sub %>% filter(metric == "PSI"),
    aes(y = PolyA_Site, fill = Time)
  ) +
    geom_col(
      aes(x = value),
      position = position_dodge(width = 0.9),
      alpha = 0.7
    ) +
    geom_errorbar(
      aes(xmin = value - sd, xmax = value + sd),
      position = position_dodge(width = 0.9),
      width = 0.2
    ) +
    geom_jitter(
      data = pts_sub %>% filter(metric == "PSI"),
      aes(x = value),
      position = position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.9
      ),
      shape = 21, color = "black", size = 1
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_fill_manual(values = custom_colors) +
    labs(x = "PSI", y = "") +
    theme_classic(base_size = 9) +
    theme(axis.text.y = element_blank())
  
  # ---------- Counts plot (no y-axis labels) ----------
  counts_plot <- ggplot(
    sum_sub %>% filter(metric == "Counts"),
    aes(y = PolyA_Site, fill = Time)
  ) +
    geom_col(
      aes(x = value),
      position = position_dodge(width = 0.9),
      alpha = 0.7
    ) +
    geom_errorbar(
      aes(xmin = value - sd, xmax = value + sd),
      position = position_dodge(width = 0.9),
      width = 0.2
    ) +
    geom_jitter(
      data = pts_sub %>% filter(metric == "Counts"),
      aes(x = value),
      position = position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.9
      ),
      shape = 21, color = "black", size = 1
    ) +
    scale_fill_manual(values = custom_colors) +
    labs(x = "Counts", y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.y = element_blank() )
  
  # ---------- combine ----------
  combined <- (psi_plot + counts_plot) + 
    plot_layout(guides = "collect") & 
    theme(
      legend.position = "none" )
  
  return(combined)
  
}

tsc2_plot <- plot_event(plotting_events[1])
tsc2_plot
syncrip_plot <- plot_event(plotting_events[2])
syncrip_plot

# PLOT TOGETHER -----------------------------------------------------------
library(cowplot)



full <- ggdraw() +
  draw_plot(overlap_type_bar,   x = 0.5, y = 0.775, width = 0.5, height = 0.225) +  #B
  
  draw_plot(genomic_distance,   x = 0, y = 0.55, width = 0.5, height = 0.225) + #C
  draw_plot(num_exons,          x = 0.5, y = 0.55, width = 0.5, height = 0.225) + #D
  
  draw_plot(gviz_grob,          x = 0, y = 0.35, width = 0.5, height = 0.18) + #E
  draw_plot(gviz_grob_SYNCRIP,    x = 0, y = 0, width = 0.5, height = 0.30) + #E
  draw_plot(tsc2_plot,           x = 0.5, y = 0.3, width = 0.5, height = 0.21) +  #F
  draw_plot(syncrip_plot,           x = 0.5, y = -0.015, width = 0.5, height = 0.3) +  #F
  
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"),
                  size  = c(12, 12, 12, 12, 12, 12),
                  x     = c(0, 0.5, 0, 0.5, 0, 0),
                  y     = c(1, 1, 0.775, 0.775, 0.55, 0.3))
full


svg(paste0("./figures/Figure_5_TSS.svg"), height = 7.5, width = 6.5)
print(full)
dev.off()


# ASD GENE ENRICHMENT -----------------------------------------
filtered_analysis_data_wide <- filtered_analysis_data_wide %>% 
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
genes <- unique(filtered_analysis_data_wide$gene_id) #1966

sig_events_df <- filtered_analysis_data_wide %>% filter(event_id %in% all_sig_events) %>% distinct(event_id, gene_id)
genes_sig <- unique(sig_events_df$gene_id) #1277

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
  
  return(c(overlap = overlap, inAnotB = inAnotB, inBnotA = inBnotA, inNeither = inNeither, p_value = fisher_result$p.value, odds_ratio = fisher_result$estimate))
}

#AS term enrichment:
#Bg genes

SFARI <- read.csv("./data/SFARI-Gene_genes_04-03-2025release_04-15-2025export.csv")

#Testing for enrichment in the each functional consequence:
SFARI.genes <- SFARI %>% filter(gene.symbol %in% genes) %>% pull(gene.symbol) %>% unique()

setA <- SFARI.genes
setB <- genes_sig
result <- compute_fisher_enrich(setA, setB, genes)
result <- as.data.frame(result)



# EXPORT SUPP. TABLES -----------------------------------------------------
both_results <- read.csv("./code/AS_APA/output/all_events_quasi_TSS.csv")
both_results <- both_results[,-c(10)]
both_emmeans <- read.csv("./code/AS_APA/output/all_events_emmeans_quasi_TSS.csv")
both_emmeans <- both_emmeans[,-c(10)]
results_list <- read.csv("./code/AS_APA/output/Fisher_result_TSS.csv")


library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "Quasi - Results", gridLines = T)
writeDataTable(wb = wb, sheet = 1, x = both_results, withFilter = F, tableStyle = "None")

addWorksheet(wb = wb, sheetName = "Quasi - EMMs", gridLines = T)
writeData(wb = wb, sheet = 2, x = both_emmeans)

addWorksheet(wb = wb, sheetName = "Fishers - Results", gridLines = T)
writeData(wb = wb, sheet = 3, x = results_list)

saveWorkbook(wb, "./tables/ATSS_AS.xlsx", overwrite = TRUE)





# PLOT HISTOGRAM OF MAX dPSI - TSS ------------------------------------------
#1) Calculate max dPSI for all sig. global events
all_sig_global_events <- both_results_all_global %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)

analysis_data_wide_summary_global_filter
max_global <- analysis_data_wide_summary_global_filter %>%
  filter(event_id %in% all_sig_global_events) %>% distinct(event_id, .keep_all = T)

global_hist <- ggplot(max_global, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Global",
       x = "<span style='font-size:9pt'>Max ΔPSI</span>",
       y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())
  
global_hist


#2) Calculate max dPSI for all sig. int. events (display both filters)
all_sig_int_events <- both_results_all_int %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)

#dPSI for each tp, across TSS sites
analysis_data_wide_summary_int_filter 
max_int1 <- analysis_data_wide_summary_int_filter %>%
  filter(event_id %in% all_sig_int_events) %>% 
  group_by(event_id) %>%
  slice_max(dpsi) %>%
  distinct(event_id, .keep_all = T) %>%
  ungroup()

max_int1_hist <- ggplot(max_int1, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Interaction",
     #  x = expression(atop("Max " * Delta * "PSI", "(per timepoint, across TSSs)")),
     x = "<span style='font-size:9pt'>Max ΔPSI</span><br>
       <span style='font-size:8pt'>(per timepoint, across TSSs)</span>",
       y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())

max_int1_hist


#dPSI for each TSS site, across tp
analysis_data_wide_summary_int_total_filter 

max_int2 <- analysis_data_wide_summary_int_total_filter %>%
  filter(event_id %in% all_sig_int_events) %>% 
  distinct(event_id, .keep_all = T)

max_int2_hist <- ggplot(max_int2, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Interaction",
      # x = expression("Max " * Delta * "PSI (per TSS, across timepoint)"),
      x = "<span style='font-size:9pt'>Max ΔPSI</span><br>
       <span style='font-size:8pt'>(per TSS, across timepoint)</span>",
      y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())

max_int2_hist

#save these files:
write.csv(max_global, "./code/AS_APA/output/max_global_plotting_TSS.csv", row.names = F)
write.csv(max_int1, "./code/AS_APA/output/max_int1_plotting_TSS.csv", row.names = F)
write.csv(max_int2, "./code/AS_APA/output/max_int2_plotting_TSS.csv", row.names = F)


##Plotting
#5'-AS

#As-3'
max_global_3 <- read.csv("./code/AS_APA/output/max_global_plotting_pA.csv")
max_int1_3 <- read.csv("./code/AS_APA/output/max_int1_plotting_pA.csv")
max_int2_3 <- read.csv("./code/AS_APA/output/max_int2_plotting_pA.csv")

global_hist_3 <- ggplot(max_global_3, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Global",
      # x = expression("Max " * Delta * "PSI"), 
       x = "<span style='font-size:9pt'>Max ΔPSI</span>",
       y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())

global_hist_3


max_int1_hist_3 <- ggplot(max_int1_3, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Interaction",
       x = "<span style='font-size:9pt'>Max ΔPSI</span><br>
       <span style='font-size:8pt'>(per timepoint, across TSSs)</span>",
       y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())


max_int1_hist_3


max_int2_hist_3 <- ggplot(max_int2_3, aes(dpsi)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Interaction",
       x = "<span style='font-size:9pt'>Max ΔPSI</span><br>
       <span style='font-size:8pt'>(per TSS, across timepoint)</span>",
       y = "# of AS Events") +
  theme_classic(base_size = 9) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x = ggtext::element_markdown())

max_int2_hist_3





#Plot together:

supp <- ggdraw() +
  draw_plot(global_hist,   x = 0, y = 0.681, width = 0.32, height = 0.319) +
  draw_plot(max_int1_hist,   x = 0.32, y = 2/3, width = 0.32, height = 1/3) +
  draw_plot(max_int2_hist,   x = 0.64, y = 2/3, width = 0.32, height = 1/3) +
  
  draw_plot(global_hist_3,   x = 0, y = 0.3477, width = 0.32, height = 0.319) +
  draw_plot(max_int1_hist_3,   x = 0.32, y = 1/3, width = 0.32, height = 1/3) +
  draw_plot(max_int2_hist_3,   x = 0.64, y = 1/3, width = 0.32, height = 1/3) +
  
  draw_plot(genes_venn,   x = 0.16, y = 0, width = 1/3, height = 1/3) +
  draw_plot(events_venn,   x = 0.57, y = 0, width = 1/3, height = 1/3) +

  draw_plot_label(label = c("A", "B", "C"),
                size  = c(12, 12, 12),
                x     = c(0, 0, 0),
                y     = c(1, 2/3, 1/3))
supp


svg(paste0("./figures/coord_supp.svg"), height = 6, width = 6)
print(supp)
dev.off()


