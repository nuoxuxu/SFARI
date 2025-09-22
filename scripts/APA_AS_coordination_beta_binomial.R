# SCRIPT: APA Alternative Splicing Analysis
# AUTHOR: Shreejoy / Gemini
# DATE: 2023-10-27
#
# DESCRIPTION:
# This script analyzes the alternative splicing of the AGO1 gene in relation to
# different polyadenylation (PolyA) sites and cell differentiation timepoints
# (iPSC, NPC, CN). It uses a beta-binomial regression model to assess the
# log-odds of exon inclusion and visualizes both the raw Percent Spliced In (PSI)
# values and the model-derived log-odds ratios.

#==============================================================================
# 1. SETUP: Load Libraries
#==============================================================================
# --- Data Manipulation and Plotting ---
library(readr)       # For reading CSV files efficiently
library(dplyr)       # For data manipulation (filter, mutate, etc.)
library(stringr)     # For string manipulation (str_sort)
library(tibble)      # For modern data frame manipulation (rownames_to_column)
library(tidyr)       # For reshaping data (pivot_wider)
library(forcats)     # For easy manipulation of factor levels (fct_rev)
library(ggplot2)     # For creating visualizations
library(cowplot)     # For arranging multiple plots into a grid
library(patchwork)   # For intuitive plot composition

# --- Modeling ---
library(glmmTMB)     # For fitting generalized linear mixed models, including beta-binomial
library(emmeans)     # For calculating estimated marginal means (model-based predictions)


#==============================================================================
# 2. DATA LOADING AND PRE-PROCESSING
#==============================================================================

# --- Load Raw Data ---
# Load the data prepared in a previous step (e.g., Python).
# The tryCatch block provides a user-friendly error if the file is not found.
tryCatch({
  raw_data <- read_csv("~/Downloads/AGO1.csv", show_col_types = FALSE)
}, error = function(e) {
  stop("Could not read the AGO1 data file. Please ensure the file path is correct.")
})

# --- Define Explicit Reference Levels for Factors ---
# Setting explicit factor levels is crucial for correct model interpretation.
# The first level of a factor is treated as the reference category in regression models.

# Order PolyA sites numerically so the most proximal site (pA1) is the reference.
pa_levels_sorted <- unique(raw_data$pas_id) %>%
  str_sort(numeric = TRUE)

# Order timepoints chronologically to treat 'iPSC' as the baseline/reference.
time_levels_sorted <- c('iPSC', 'NPC', 'CN')

#' Filters polyA sites based on their proportional usage across the entire dataset.
#' This is a crucial step to remove noise from infrequently used sites.
#'
#' @param data_frame The raw, long-format data frame from the input CSV.
#' @param threshold Numeric (0-1). pA sites with total usage below this are removed.
#' @return A filtered data frame.
filter_pa_sites_by_usage <- function(data_frame, threshold = 0.025) {
  # Calculate total reads per pA site and the grand total
  pa_site_totals <- data_frame %>%
    group_by(pas_id) %>%
    summarise(site_total_reads = sum(total), .groups = 'drop')
  
  grand_total_reads <- sum(pa_site_totals$site_total_reads)
  
  # Identify which pA sites meet the usage threshold
  sites_to_keep <- pa_site_totals %>%
    mutate(proportion = site_total_reads / grand_total_reads) %>%
    filter(proportion >= threshold) %>%
    pull(pas_id) # Extract names of sites to keep
    
  cat("Identified", length(sites_to_keep), "pA sites meeting the >=", threshold*100, "% usage threshold:\n")
  print(sites_to_keep)
  
  # Return the original data frame filtered for the selected sites
  return(data_frame %>% filter(pas_id %in% sites_to_keep))
}

#' Filters polyA sites based on a minimum read count across a minimum number of replicates.
#'
#' This is a robust filtering method to ensure that a site is sufficiently
#' expressed across replicates to be included in statistical modeling.
#'
#' @param data_frame The raw, long-format data frame from the input CSV.
#'                  It must contain columns for 'pas_id', 'timepoint', and 'total'.
#' @param min_count The minimum number of reads required in a sample for it to count.
#' @param min_replicates The minimum number of biological replicates that must meet the min_count.
#'
#' @return A filtered data frame.

filter_pa_sites_by_replicates <- function(data_frame, min_count = 5, min_replicates = 3) {
  
  # 1. For each pA site, count how many replicates meet the minimum read threshold.
  sites_to_keep <- data_frame %>%
    # First, identify which samples meet the read count threshold for each site
    filter(total >= min_count) %>%
    # Now, group by pA site and count the number of unique replicates
    group_by(pas_id) %>%
    summarise(n_replicates_meeting_threshold = n_distinct(timepoint), .groups = 'drop') %>%
    # Finally, keep only the pA sites that meet the minimum replicate threshold
    filter(n_replicates_meeting_threshold >= min_replicates) %>%
    pull(pas_id) # Extract the names of the sites to keep
    
  # 2. Filter the original data frame and return it
  cat("Identified", length(sites_to_keep), "pA sites with at least", min_count, "reads in at least", min_replicates, "replicates:\n")
  print(sites_to_keep)
  
  return(data_frame %>% filter(pas_id %in% sites_to_keep))
}

# --- Filter Data for Robustness ---
# This custom function (assumed to be defined elsewhere) filters out PolyA sites
# that do not have a minimum read count across a minimum number of replicates.
# This step increases confidence that the analyzed sites are reliably expressed.
# NOTE: The specific thresholds (min_count, min_replicates) are user-defined.
filtered_raw_data <- filter_pa_sites_by_replicates(raw_data, min_count = 50, min_replicates = 5)

# --- Prepare Analysis Data Frame ---
# Apply the defined factor levels and prepare the main data frame for analysis.
analysis_data <- filtered_raw_data %>%
  mutate(
    BioSampleID = factor(timepoint),
    PolyA_Site = factor(pas_id, levels = pa_levels_sorted),
    Time = factor(sub("_.*", "", BioSampleID), levels = time_levels_sorted)
  )

# --- Reshape Data to Wide Format ---
# The beta-binomial model in glmmTMB requires a two-column response variable:
# one for successes (Included reads) and one for failures (Skipped reads).
# We pivot the data from a long format to this wide format.
analysis_data_wide <- analysis_data %>%
  pivot_wider(
    names_from = type,      # Create new columns named "Included" and "Skipped"
    values_from = total,    # Populate these columns with values from the 'total' column
    values_fill = 0         # Fill missing combinations with 0
  )

# --- Add Pseudocount and Calculate PSI ---
# A pseudocount of 1 is added to both Included and Skipped reads.
# This prevents division-by-zero errors when calculating PSI and helps stabilize
# variance in samples with very low counts, which improves model stability.
analysis_data_wide <- analysis_data_wide %>%
  mutate(
    Included = Included + 1,
    Skipped = Skipped + 1,
    psi = Included / (Included + Skipped) # Calculate Percent Spliced In (PSI) for each sample
  )

#==============================================================================
# 3. BETA-BINOMIAL MIXED-EFFECTS MODELING AND HYPOTHESIS TESTING
#==============================================================================

# --- Fit Nested Generalized Linear Mixed-Effects Models (GLMMs) ---
# We fit three nested models to test our hypotheses using Likelihood Ratio Tests (LRT).
# A random intercept `(1 | BioSampleID)` is included in all models to account for
# the non-independence of measuring multiple PolyA sites from the same biological sample.

# Model 1: Interaction Model (Time * PolyA_Site)
# Tests if the effect of PolyA site on splicing differs across timepoints,
# while controlling for sample-to-sample variability.
bb_fit_interaction <- glmmTMB(
  cbind(Included, Skipped) ~ Time * PolyA_Site + (1 | BioSampleID),
  data = analysis_data_wide,
  family = betabinomial()
)

# Model 2: Additive Model (Time + PolyA_Site)
# Assumes the effects of Time and PolyA site are independent, while controlling
# for sample-to-sample variability.
bb_fit_full <- glmmTMB(
  cbind(Included, Skipped) ~ Time + PolyA_Site + (1 | BioSampleID),
  data = analysis_data_wide,
  family = betabinomial()
)

# Model 3: Reduced Model (Time only)
# Serves as a baseline, accounting only for the effect of Time and
# sample-to-sample variability.
bb_fit_reduced <- glmmTMB(
  cbind(Included, Skipped) ~ Time + (1 | BioSampleID),
  data = analysis_data_wide,
  family = betabinomial()
)

# --- Perform Likelihood Ratio Tests (LRT) ---
# The LRT compares nested models to determine if the additional fixed-effect terms
# provide a statistically significant improvement in fit, given the random effect structure.

# Test 1: Global effect of PolyA site.
# Compares the additive model to the reduced model. A significant p-value
# suggests that PolyA_Site as a whole is a significant predictor of splicing.
lrt_result_global <- anova(bb_fit_reduced, bb_fit_full)
print("Likelihood Ratio Test for the Global polyA Term:")
print(lrt_result_global)

# Test 2: Significance of the interaction term.
# Compares the interaction model to the additive model. A significant p-value
# indicates that the relationship between PolyA site and splicing is different
# at different timepoints.
lrt_result_interaction <- anova(bb_fit_interaction, bb_fit_full)
print("Likelihood Ratio Test for the Interaction polyA*Time Term:")
print(lrt_result_interaction)

# --- Extract and Save Model Coefficients ---
# Extract the fixed-effects coefficients (beta coefficients) from the most
# complex model to examine the estimated effects of each predictor.
coef_table <- summary(bb_fit_interaction)$coefficients$cond

# Convert the matrix of coefficients to a tidy data frame.
# The 'Term' column will contain the predictor names (e.g., TimeNPC, PolyA_SitepA2).
# The 'Estimate' column contains the log-odds ratio for that term relative to the reference.
coef_df <- as.data.frame(coef_table) %>%
  rownames_to_column(var = "Term")

# Save the coefficients to a CSV file for supplementary materials or further analysis.
write_csv(coef_df, "AGO1_interaction_model_coefficients.csv")

# Print a confirmation message to the console.
print("Model coefficients have been saved to AGO1_interaction_model_coefficients.csv")


#==============================================================================
# 4. VISUALIZATION
#==============================================================================

# --- Define Custom Color Palette ---
# Define colors for each timepoint to ensure consistency across plots.
# These hex codes are selected to match the original publication figure.
custom_colors <- c(
  "iPSC" = "#C1B2E0",
  "NPC"  = "#907FB9",
  "CN"   = "#3D3B8E"
)


# --- Plot 1: Percent Spliced In (PSI) ---

# Prepare data for plotting. Factor levels are reversed (`fct_rev`) so that
# they appear in the desired order on the y-axis (e.g., pA1, pA2, ... from top to bottom).
plot_data <- analysis_data_wide %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("iPSC", "NPC", "CN"))),
    PolyA_Site = fct_rev(as.factor(PolyA_Site))
  )

# Calculate mean PSI for each group to determine the height of the bars.
summary_data <- plot_data %>%
  group_by(PolyA_Site, Time) %>%
  summarise(mean_psi = mean(psi, na.rm = TRUE), .groups = 'drop')

# Generate the PSI plot.
psi_plot <- ggplot(summary_data, aes(y = PolyA_Site, fill = Time)) +
  # Bars represent the mean PSI value for each group.
  geom_col(
    aes(x = mean_psi),
    position = position_dodge(width = 0.9),
    alpha = 0.7 # Use transparency to make overlaid points more visible.
  ) +
  # Points represent the PSI value for each individual biological replicate.
  geom_jitter(
    data = plot_data,
    aes(x = psi),
    position = position_jitterdodge(
      jitter.width = 0.2, # Control horizontal spread of points.
      dodge.width = 0.9   # Ensure points align with their corresponding bars.
    ),
    shape = 21, color = "black", size = 2
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "PSI (Percent Spliced In)",
    y = "PolyA Site"
  ) +
  theme_cowplot() +
  theme(legend.position = "none") # Hide legend for this plot; it will be shown on the combined plot.


# --- Plot 2: Log Odds Ratio from Model ---

# Use emmeans to calculate the estimated marginal means from the interaction model.
# This is the standard way to get model-based estimates (log-odds) and confidence
# intervals for each group, correctly accounting for all model terms.
emm_fit <- emmeans(bb_fit_interaction, specs = ~ Time | PolyA_Site)
emm_df <- as.data.frame(emm_fit)

# Prepare the emmeans data for plotting, reversing factor levels as before.
emm_df_ordered <- emm_df %>%
  mutate(
    Time = fct_rev(factor(Time, levels = c("iPSC", "NPC", "CN"))),
    PolyA_Site = fct_rev(as.factor(PolyA_Site))
  )

# Generate the Log Odds Ratio plot.
log_or_plot <- ggplot(emm_df_ordered, aes(x = emmean, y = PolyA_Site, fill = Time)) +
  # Bars represent the estimated log-odds ratio from the model.
  geom_col(
    position = position_dodge(width = 0.9),
    alpha = 0.7 # Match transparency of the PSI plot.
  ) +
  # Error bars represent the 95% confidence interval for the estimate.
  geom_errorbar(
    aes(xmin = asymp.LCL, xmax = asymp.UCL),
    position = position_dodge(width = 0.9),
    width = 0.25
  ) +
  # A vertical dashed red line at x=0 indicates no change in odds.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Log Odds Ratio (Inclusion vs. Exclusion) \u00B1 95% CI", # \u00B1 is the +/- symbol
    y = NULL, # Remove y-axis label to avoid redundancy in the final combined plot.
    fill = "Cell Stage"
  ) +
  theme_cowplot() +
  guides(fill = guide_legend(reverse = TRUE)) +

  # Remove y-axis text and ticks for a cleaner look when combined.
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#==============================================================================
# 5. COMBINE AND DISPLAY PLOTS
#==============================================================================

# Arrange the PSI plot and the Log Odds Ratio plot side-by-side into a single figure.

# --- Construct a Dynamic Title ---
# Extract the event ID, assuming it's the same for all rows in the dataset.
event_identifier <- analysis_data_wide$event_id[1]

# Extract p-values from the LRT results (index [2] corresponds to the p-value row).
p_global <- lrt_result_global$`Pr(>Chisq)`[2]
p_interaction <- lrt_result_interaction$`Pr(>Chisq)`[2]

# Format the p-values for clean display in the title.
p_global_formatted <- format.pval(p_global)
p_interaction_formatted <- format.pval(p_interaction)

# Create the final title string. Using a newline character `\n` for better spacing.
plot_main_title <- paste0(
  event_identifier,
  "\nGlobal p = ", p_global_formatted,
  ", Interaction p = ", p_interaction_formatted
)


# --- Arrange Plots and Add Title ---
# Arrange the plots side-by-side using patchwork.
final_plot <- psi_plot + log_or_plot +
  plot_layout(widths = c(0.8, 1)) + # Give slightly more width to the plot with the legend.
  plot_annotation(
    title = plot_main_title,
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12)) # Center the title
  )

# Display the final combined plot.
print(final_plot)
