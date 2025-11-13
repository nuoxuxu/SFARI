# DESCRIPTION:
# This script analyzes the alternative splicing in relation to
# different polyadenylation (PolyA) sites and cell differentiation timepoints
# (t00, t04 and t30). It uses a quasi-binomial regression model to assess the
# log-odds of exon inclusion and visualizes the results.

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



setwd("")

# --- Load Raw Data ---
raw_data <- read_csv("./code/AS_APA/input/pA_counts.csv", show_col_types = FALSE)




#Function to filter polyA sites
filter_pa_sites_by_replicates_new <- function(data_frame, min_count = 5, min_replicates = 3, min_count_all_reps = 1) {
  
  # Step 1: Filter pas_ids that meet min_count_all_reps in *every* sample (ignoring type)
  if (!is.null(min_count_all_reps)) {
    pas_ids_to_keep <- data_frame %>%
      group_by(event_id, pas_id, timepoint) %>%
      summarise(total_sum = sum(total, na.rm = TRUE), .groups = "drop") %>%
      group_by(event_id, pas_id) %>%
      summarise(all_timepoints_ok = all(total_sum >= min_count_all_reps), .groups = "drop") %>%
      filter(all_timepoints_ok) %>%
      dplyr::select(event_id, pas_id)
    
    data_frame <- data_frame %>%
      semi_join(pas_ids_to_keep, by = c("event_id", "pas_id"))
    
    cat("Retained", nrow(pas_ids_to_keep), "pA sites passing min_count_all_reps =", min_count_all_reps, "\n")
  }
  # Step 2: Filter rows that meet the read count threshold
  filtered <- data_frame %>%
    filter(total >= min_count)
  
  # Step 3: For each event_id and pas_id, count how many unique timepoints meet the threshold
  sites_to_keep <- filtered %>%
    group_by(event_id, pas_id) %>%
    summarise(
      n_replicates = n_distinct(timepoint),
      .groups = 'drop'
    ) %>%
    filter(n_replicates >= min_replicates)
  
  cat("Filtered to", nrow(sites_to_keep), "pA sites across", 
      n_distinct(sites_to_keep$event_id), "event IDs meeting thresholds.\n")
  
  # Step 4: Use semi_join to keep only the matching rows from the original data
  filtered_data <- data_frame %>%
    semi_join(sites_to_keep, by = c("event_id", "pas_id"))
  
  return(filtered_data)
}





filtered_raw_data <- filter_pa_sites_by_replicates_new(raw_data, 
                                                       min_count = 50, 
                                                       min_replicates = 5, 
                                                       min_count_all_reps = 10) 
# Retained 52470 pA sites passing min_count_all_reps = 10 
# Filtered to 45918 pA sites across 22317 event IDs meeting thresholds.

# --- Prepare Analysis Data Frame ---
pa_levels_sorted <- unique(raw_data$pas_id) %>%
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
    PolyA_Site = factor(pas_id, levels = pa_levels_sorted),
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


write.csv(analysis_data_wide, "./code/AS_APA/output/pA_counts_wide_new_filter.csv", 
          row.names = F)





# PREPARE EVENTS ----------------------------------------------------------
analysis_data_wide <- read.csv("./code/AS_APA/output/pA_counts_wide_new_filter.csv")
#Filter to event_ids w/ >1 pA
filtered_analysis_data_wide <- analysis_data_wide %>%
  group_by(event_id) %>%
  filter(n_distinct(pas_id) > 1) %>%
  ungroup()
length(unique(filtered_analysis_data_wide$event_id)) #11,766
  
#Categorize events to 3 groups (those that pass global only, pass interaction only, and both)

#1) Find events that pass global dPSI:
analysis_data_wide_summary_global <- filtered_analysis_data_wide %>%
  group_by(event_id, pas_id) %>%
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


#2) Find events that pass interaction dPSI: (test delta for each tp, across pA site)
analysis_data_wide_summary_int <- filtered_analysis_data_wide %>%
  group_by(event_id, pas_id, Time) %>%
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
  group_by(event_id, pas_id, Time) %>%
  summarise(mean_psi = mean(psi)) %>% 
  ungroup() %>%
  group_by(event_id, pas_id) %>%
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
global_only_events <- setdiff(filtered_analysis_data_wide_global$event_id, filtered_analysis_data_wide_int_v2$event_id) #991
int_only_events <- setdiff(filtered_analysis_data_wide_int_v2$event_id, filtered_analysis_data_wide_global$event_id) #319
both_events <- intersect(filtered_analysis_data_wide_int_v2$event_id, filtered_analysis_data_wide_global$event_id) #1677

all_global <- unique(filtered_analysis_data_wide_global$event_id) #2668
all_int <- unique(filtered_analysis_data_wide_int_v2$event_id) #1996



# QUASI BINOMIAL -----------------------------------------------------
events <- unique(filtered_analysis_data_wide$event_id)
test <- filtered_analysis_data_wide %>% filter(event_id %in% events)


fit_both_models_and_extract <- function(event_data, event_id) {
  tryCatch({
    # --- Preprocessing ---
    event_data <- event_data %>%
      mutate(
        Time = factor(Time, levels = c("t00", "t04", "t30")),
        PolyA_Site = factor(PolyA_Site, levels = str_sort(unique(PolyA_Site), numeric = TRUE))
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
        cbind(Included, Skipped) ~ Time + PolyA_Site,
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
        cbind(Included, Skipped) ~ Time * PolyA_Site,
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
    emm_fit <- emmeans(fit_interaction, specs = ~ Time | PolyA_Site)
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
      PolyA_Site = NA,
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
length(unique(results_pass$event_id)) #11750
write.csv(results_pass, "./code/AS_APA/output/all_events_quasi.csv", row.names = F)


emmeans_df <- emmeans_df[,-c(9,10)]
emmeans_pass <- emmeans_df %>% filter(interaction_model_converged == T)
length(unique(emmeans_pass$event_id)) #11750
write.csv(emmeans_pass, "./code/AS_APA/output/all_events_emmeans_quasi.csv", row.names = F)



# ANALYZING EVENTS --------------------------------------------------------
# Apply effect size filter:
both_results <- read.csv("./code/AS_APA/output/all_events_quasi.csv")

both_results %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #9826
both_results %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #3835
both_results %>% filter(global_LRT_FDR <= 0.05 &
                          interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #3463

both_results_filter <- both_results %>% filter(event_id %in% both_events) #1670
both_results_filter %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #1668
both_results_filter %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #647
both_results_filter %>% filter(global_LRT_FDR <= 0.05 &
                          interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #646


both_results_all_global <- both_results %>% filter(event_id %in% all_global) 
both_results_all_global %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #2650


both_results_all_int <- both_results %>% filter(event_id %in% all_int)
both_results_all_int %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% tally() #785


# Compare sig events to sig events from Fisher's - can only compare global
results_list <- read.csv("./code/AS_APA/output/Fisher_result.csv")

results_list <- results_list %>%
  mutate(
    gene_id = sub(";.*", "", event_id)
    
  )
global_results <- both_results %>% filter(event_id %in% all_global) 
global_results_events <- global_results %>% filter(global_LRT_FDR <= 0.05)
  length(unique(global_results_events$event_id)) #2650 sig global events
global_results_events_unique <- unique(global_results_events$event_id)
#   
#Here are the Fisher's events w/ at least one padj < 0.05 and |logOR|>0.05
sig_results #2813 unique
sig_results_unique <- unique(sig_results$event_id)

##But the 'universe' of testable events is different in each:
#Events passing read counts filters for quasi:
quasi_universe <- unique(filtered_analysis_data_wide$event_id) #11766
fisher_universe <- unique(results_list$event_id) #3503

overlapping_universe <- intersect(quasi_universe, fisher_universe) #3169

global_results_events_unique2 <- global_results_events_unique[global_results_events_unique %in% overlapping_universe] #1771
sig_results_unique2 <- sig_results_unique[sig_results_unique %in% overlapping_universe] #2569

length(intersect(global_results_events_unique2, sig_results_unique2)) #1704
length(setdiff(global_results_events_unique2, sig_results_unique2)) #67
length(setdiff(sig_results_unique2, global_results_events_unique2)) #865

#Thus, of the events tested in both analyses (i.e. passing read count filters) = 3169 total
  #Then, filtered to sig events
  #1704 overlapping (out of 1704+67+865 = 2636 total sig)
  #65%




# PLOTTING ----------------------------------------------------
# Load in data
both_results <- read.csv("./code/AS_APA/output/all_events_quasi.csv")
both_emmeans <- read.csv("./code/AS_APA/output/all_events_emmeans_quasi.csv")

filtered_df_ES <- read.csv("./code/AS_APA/input/transcript_pA_sites.csv")
colnames(filtered_df_ES)[5] <- "pas_id"

counts_annotated <- read.csv("./code/AS_APA/input/tr_pA_counts.csv")


#All events with at least 1 type of coordination: #unique set of global + int
both_results_all_global <- both_results %>% filter(event_id %in% all_global) 
all_sig_global_events <- both_results_all_global %>% filter(global_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)


both_results_all_int <- both_results %>% filter(event_id %in% all_int)
all_sig_int_events <- both_results_all_int %>% filter(interaction_LRT_FDR <= 0.05) %>% distinct(event_id) %>% pull(event_id)

all_sig_events <- unique(union(all_sig_global_events, all_sig_int_events) ) #2789


#Final event_table:
event_table <- filtered_analysis_data_wide %>% filter(event_id %in% all_sig_events)

# GENOMIC DISTANCE --------------------------------------------------------
#To get proximal coord per event_id-pA site:
unique_prox <- filtered_df_ES %>% distinct(gene_id, pas_id, cluster_proximal_coord)

#To get all event_ids to polyA sites:
event_table <- event_table %>% distinct(event_id, pas_id) %>%
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
event_table <- merge(event_table, unique_prox, 
                     by.x = c("gene_id", "pas_id"), by.y = c("gene_id", "pas_id"))


#Extract exon coords:
coords <- event_table %>%
   # Extract only the coordinate part of event_id: skip everything up to the first colon (after chr...)
  mutate(coord_part = str_replace(event_id, "^[^:]+:[^:]+:", ""),  # removes up to second colon
         
         # Now extract just the genomic numbers (this avoids chromosome number)
         coord_nums = str_extract_all(coord_part, "\\d+"),
         coord_nums = lapply(coord_nums, as.numeric), 
         strand = str_extract(event_id, ".$") ) %>%
  
  # Calculate most 3' coord based on strand
  rowwise() %>%
  mutate(
    most_3prime_coord = if (strand == "+") max(coord_nums, na.rm = TRUE) else min(coord_nums, na.rm = TRUE),
    distance_to_3prime = abs(cluster_proximal_coord - most_3prime_coord)
  ) %>%
  ungroup()

coords <- coords %>% distinct(event_id, pas_id, .keep_all = T)


#Filter to the most proximal pA per event_id:
coords_filtered <- coords %>% 
  mutate(pA_num = as.integer(sub("pA", "", pas_id))   # extract numeric part
  ) %>%
  group_by(event_id) %>%
  slice_min(order_by = pA_num, with_ties = FALSE) %>%
  ungroup() #2789

median_distance <- median(coords_filtered$distance_to_3prime) #12,527


genomic_distance <- coords_filtered %>%
  ggplot(aes(x = distance_to_3prime )) +
  geom_density(alpha = 0.35, linewidth = 0.6, trim = T, fill = "cornflowerblue", colour = "darkblue") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(colour = "black"), 
        axis.title.x = element_markdown() ) +
  scale_x_log10(
    labels = trans_format("log10", math_format(10^.x)), 
    limits = c(9, NA)
  ) +
  geom_vline(xintercept = median_distance, colour = "black", linetype = "dashed") +
  labs(
    x = "<span style='font-size:9pt'>Genomic Distance (nt)</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
    y = "Density")

genomic_distance
 



# NUMBER OF EXONS ---------------------------------------------------------
ORFanage_replaced <- rtracklayer::import("./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_gtf.gtf") #only using 158,844 here

#Filter to only transcripts in the significant event-pA IDs (only for the most proximal pA)
counts_filter <- counts_annotated %>% 
  filter(paste(event_id, pas_id) %in% paste(coords_filtered$event_id, coords_filtered$pas_id)) #2813
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
  ungroup() #2789


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
         selected_coord = map2_dbl(coord_nums, strand, ~ if(.y == "+") .x[3] else .x[2]) )


#Now get # exons b/w:
exons <- ORFanage_replaced[ORFanage_replaced$type == "exon"]

# Step 1: create regions GRanges for each event_id
# Make sure start <= end even for negative strand
regions <- GRanges(
  seqnames = as.character(counts_max$chr),
  ranges   = IRanges(
    start = pmin(counts_max$selected_coord, counts_max$transcript_end),
    end   = pmax(counts_max$selected_coord, counts_max$transcript_end)
  ),
  strand = as.character(counts_max$strand),
  transcript_id = counts_max$transcript_id,
  event_id      = counts_max$event_id
)

# Step 2: count exons per region (event_id)
#Include exons that overlap transcript_end (right boundary).
#Exclude exons that overlap selected_coord (left boundary).

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
        start = start(region_i) + 1, 
        end   = end(region_i)
      ),
      strand   = strand(region_i)
    )
  } else {
    region_i_adj <- GRanges(
      seqnames = seqnames(region_i),
      ranges   = IRanges(
        start = start(region_i), 
        end   = end(region_i) - 1
      ),
      strand   = strand(region_i)
    )
  }
  
  # count overlaps
  length(findOverlaps(region_i_adj, exons_tx, ignore.strand = FALSE))
})

# Step 3: assign back to counts_max
counts_max$num_exons_in_region <- num_exons_in_region

median_num_exon <- median(counts_max$num_exons_in_region) #5

num_exons <- counts_max %>%  
  ggplot(aes(x = num_exons_in_region)) +
  geom_bar(width = 0.7, alpha = 0.35, linewidth = 0.25,
           fill = "cornflowerblue", color = "darkblue") +
  theme_classic(base_size = 9) +
  coord_cartesian(ylim = c(0, 375)) +
  labs(
    x = "<span style='font-size:9pt'>Number of Exons</span><br><span style='font-size:8pt'>(as depicted in A)</span>",
    y = "# of AS events") +
  theme(axis.text.x = element_text(color = "black"), 
        axis.title.x = element_markdown(),
        axis.text.y = element_text(colour = "black"), 
        axis.ticks = element_line(colour = "black") ) +
  geom_vline(xintercept = median_num_exon, colour = "black", linetype = "dashed") 

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
  ungroup() #2789


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
gtf_sub <- ORFanage_replaced[ORFanage_replaced$transcript_id %in% unique(max_tr$transcript_id)] #should be only 2415 transcripts

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



# PLOTTING AGO1 -----------------------------------------------------------
AGO1 <- counts_annotated %>% filter(event_id == "AGO1;SE:chr1:35888610-35892557:35892677-35893097:+") %>%
  distinct(event_id, gene_id, transcript_id, pas_id) %>%
  filter(pas_id %in% c("pA4", "pA8", "pA10"))


#Plot 3 sep plots, 1 for each timepoint:
timepoints <- c("sum_t00", "sum_t04", "sum_t30")
pas_ids <- AGO1 %>% distinct(pas_id) %>% pull(pas_id)
desired_order <- c("pA4", "pA8", "pA10")
pas_ids <- desired_order[desired_order %in% pas_ids]



AGO1$pas_id[AGO1$pas_id == "pA4"] <- "pA1"
AGO1$pas_id[AGO1$pas_id == "pA8"] <- "pA2"
AGO1$pas_id[AGO1$pas_id == "pA10"] <- "pA3"

pas_ids <- AGO1 %>% distinct(pas_id) %>% pull(pas_id)
desired_order <- c("pA1", "pA2", "pA3")
pas_ids <- desired_order[desired_order %in% pas_ids]


AGO_gtf_309 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.1207.309") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(AGO_gtf_309, "./code/AS_APA/input/AGO_gtf_309.gtf")

AGO_gtf_317 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.1207.317") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(AGO_gtf_317, "./code/AS_APA/input/AGO_gtf_317.gtf")

AGO_gtf_336 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.1207.336") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(AGO_gtf_336, "./code/AS_APA/input/AGO_gtf_336.gtf")


txdb_309 <- makeTxDbFromGFF("./code/AS_APA/input/AGO_gtf_309.gtf", format = "gtf")
txdb_317 <- makeTxDbFromGFF("./code/AS_APA/input/AGO_gtf_317.gtf", format = "gtf")
txdb_336 <- makeTxDbFromGFF("./code/AS_APA/input/AGO_gtf_336.gtf", format = "gtf")


fill_colors <- rep("steelblue", 21) 
fill_colors[c(4)] <- "red" ##This is manual...

gviz_fontsize <- 7
# 

txTr_309 <- GeneRegionTrack(txdb_309,
                            col = NA,
                            fill = fill_colors,
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA1")
displayPars(txTr_309) <- list(transcript = "pA4",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_317 <- GeneRegionTrack(txdb_317,
                            col = NA,
                            fill = fill_colors,
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA2")
displayPars(txTr_317) <- list(transcript = "pA8",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_336 <- GeneRegionTrack(txdb_336,
                            col = NA,
                            fill = fill_colors,
                            showId = F, 
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA3")
displayPars(txTr_336) <- list(transcript = "pA10",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 



# Save Gviz as grob
gviz_grob <- grid::grid.grabExpr(
  plotTracks(c(GenomeAxisTrack(name = "", scale = 0.25), 
               txTr_309, txTr_317, txTr_336),
             # main = "AGO1 at t00",
             showTitle = T,
             cex.main = 0.8,
             just.group = "left",
             background.title = "white",
             col.title        = "black",
             background.panel = "white",
             col.axis         = "white", #hack to hide the axis
             col.frame        = "white", 
             sizes = c(0.1, 0.1, 0.1, 0.1), #0.925
             title.width = 0.5)
)



# PLOTTING SMARCA4 ---------------------------------------------------------
#Representative tr is the most highly expressed:
SMARCA4 <-  counts_annotated %>% filter(event_id == "SMARCA4;SE:chr19:11033517-11033767:11033865-11034123:+") %>%
  distinct(event_id, gene_id, transcript_id, pas_id) %>%
  filter(pas_id %in% c("pA4", "pA5", "pA6"))

#Representative tr is the most highly expressed:
# "PB.100786.150" "PB.100786.175" "PB.100786.185"

SMARCA4_gtf_150 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.100786.150") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SMARCA4_gtf_150, "./code/AS_APA/input/SMARCA4_gtf_150.gtf")

SMARCA4_gtf_175 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.100786.175") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SMARCA4_gtf_175, "./code/AS_APA/input/SMARCA4_gtf_175.gtf")

SMARCA4_gtf_185 <- ORFanage_replaced %>%
  subset(
    (mcols(ORFanage_replaced)$transcript_id == "PB.100786.185") &
      (mcols(ORFanage_replaced)$type != "start_codon")
  )
rtracklayer::export(SMARCA4_gtf_185, "./code/AS_APA/input/SMARCA4_gtf_185.gtf")

txdb_150 <- makeTxDbFromGFF("./code/AS_APA/input/SMARCA4_gtf_150.gtf", format = "gtf")
txdb_175 <- makeTxDbFromGFF("./code/AS_APA/input/SMARCA4_gtf_175.gtf", format = "gtf")
txdb_185 <- makeTxDbFromGFF("./code/AS_APA/input/SMARCA4_gtf_185.gtf", format = "gtf")


fill_colors <- rep("steelblue", 37) 
fill_colors[c(28)] <- "red" ##This is manual...

gviz_fontsize <- 7

txTr_150 <- GeneRegionTrack(txdb_150,
                            col = NA,
                            fill = "steelblue",
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA1")
displayPars(txTr_150) <- list(transcript = "pA1",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_175 <- GeneRegionTrack(txdb_175,
                            col = NA,
                            fill = fill_colors,
                            fontcolor.group = "transparent",
                            showId = F,
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA2"
)
displayPars(txTr_175) <- list(transcript = "pA2",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.height = 25*0.3) 

txTr_185 <- GeneRegionTrack(txdb_185,
                            col = NA,
                            fill = "steelblue",
                            showId = F, 
                            collapse = FALSE,
                            min.distance = 0,
                            name = "pA3"
)
displayPars(txTr_185) <- list(transcript = "pA3",
                              showId = FALSE, 
                              fontsize = gviz_fontsize, 
                              min.distance = 0,
                              min.height = 25*0.3) 



# Save Gviz as grob
gviz_grob_SMARCA4 <- grid::grid.grabExpr(
  plotTracks(c(GenomeAxisTrack(name = "", scale = 0.25), 
               txTr_150, txTr_175, txTr_185),
             # main = "SMARCA4 at t00",
             from = 11020515, 
             to = 11062666,
             showTitle = T,
             cex.main = 0.8,
             just.group = "left",
             background.title = "white",
             col.title        = "black",
             background.panel = "white",
             col.axis         = "white", #hack to hide the axis
             col.frame        = "white", 
             sizes = c(0.1, 0.1, 0.1, 0.1), #0.925
             title.width = 0.5)
)


# PLOTTING AGO1 & SMARCA4 PSI -----------------------------------------------
plotting_events <- c("AGO1;SE:chr1:35888610-35892557:35892677-35893097:+", 
                     "SMARCA4;SE:chr19:11033517-11033767:11033865-11034123:+")
pa_levels_sorted <- unique(raw_data$pas_id) %>%
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
                position = position_dodge(width = 0.9), 
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
  coord_cartesian(xlim = c(0, 1)) + 
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




# PLOT TOGETHER -----------------------------------------------------------
library(cowplot)

full <- ggdraw() +
  draw_plot(overlap_type_bar,   x = 0.5, y = 0.775, width = 0.5, height = 0.225) +  #B
  
  draw_plot(genomic_distance,   x = 0, y = 0.55, width = 0.5, height = 0.225) + #C
  draw_plot(num_exons,          x = 0.5, y = 0.55, width = 0.5, height = 0.225) + #D

  draw_plot(gviz_grob,          x = 0, y = 0.275, width = 0.61, height = 0.275) +
  draw_plot(gviz_grob_SMARCA4,    x = 0, y = 0, width = 0.61, height = 0.275) + 
  draw_plot(psi_plot,           x = 0.6, y = 0, width = 0.4, height = 0.55) +  
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"),
                  size  = c(12, 12, 12, 12, 12, 12),
                  x     = c(0, 0.5, 0, 0.5, 0, 0),
                  y     = c(1, 1, 0.775, 0.775, 0.55, 0.275))
full



svg(paste0("./figures/Figure_5.svg"), height = 7.5, width = 6.5)
print(full)
dev.off()


# ASD GENE ENRICHMENT -----------------------------------------
filtered_analysis_data_wide <- filtered_analysis_data_wide %>% 
  mutate(gene_id = str_extract(event_id, "^[^;]+"))
genes <- unique(filtered_analysis_data_wide$gene_id) #3392

sig_events_df <- filtered_analysis_data_wide %>% filter(event_id %in% all_sig_events) %>% distinct(event_id, gene_id)
genes_sig <- unique(sig_events_df$gene_id) #1606

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


