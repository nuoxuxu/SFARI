# Transcription Start Site (TSS) Coordination Analysis

This pipeline quantifies differential transcription start site (TSS) usage across a three neuronal differentiation timepooints (t00, t04 and t30).

> **Note:** The same framework is applied to alternative splicing (AS) and alternative polyadenylation (APA) coordination by substituting the appropriate count matrix and reformulating the PSI-equivalent metric. The modelling code is agnostic to event type and requires only a consistent `Included`/`Skipped` long-format input.

---

## 1. Input data

The raw input (`TSS_counts.csv`) is in long format with one row per combination of `event_id`, `TSS_id`, `timepoint`, and `type` (Included/Skipped). The `total` column holds read counts. 

---

## 2. Quality filtering

Filtering is applied in two sequential steps inside `filter_TSS_sites_by_replicates_new()`.

### Step 1 — Baseline expression filter

For each `(event_id, TSS_id)` pair, reads are summed across Included/Skipped within each timepoint. Any site with fewer than `min_count_all_reps = 10` summed reads at *any* timepoint is dropped entirely. This ensures that usage estimates are not computed for sites that are essentially absent at one timepoint. 

### Step 2 — Replicate-depth filter

Individual sample rows are then required to have `total ≥ 50`, and a site must meet this threshold in at least `min_replicates = 5` samples. This removes sites that pass the timepoint-level sum only because a single high-count replicate inflates the aggregate.

---

## 3. PSI calculation

Data are reshaped to wide format (one row per `event_id × TSS_id × sample`), with separate columns for Included and Skipped counts. A pseudocount of 1 is added to both before computing:

```
PSI = Included / (Included + Skipped)
```

The pseudocount prevents undefined PSI values. 

---

## 4. Event categorisation

Only events with more than one detected `TSS_id` are eligible for differential TSS analysis. These are stratified into three categories based on mean PSI summaries:

| Category | Definition
|---|---|---|
| **Global** | ∆PSI between TSS sites (max − min), averaged across timepoints, ≥ 0.1
| **Interaction** | ∆PSI across timepoints for at least 1 TSS site ≥ 0.1 *and* ∆PSI between TSS sites at at least 1 timepoint ≥ 0.1 
| Global only | Global threshold met, interaction threshold not met 
| Interaction only | Interaction threshold met, global threshold not met 
| Both | Both thresholds met 

---

## 5. Quasi-binomial GLM

For each event, three nested GLMs are fit with a quasi-binomial family (which accounts for overdispersion in count ratios without requiring explicit estimation of a dispersion parameter):

| Model | Formula | Purpose |
|---|---|---|
| Reduced | `~ Time` | Baseline time effect |
| Full | `~ Time + TSS` | Adds TSS main effect |
| Interaction | `~ Time * TSS` | Adds TSS × time interaction |

Two likelihood ratio tests (F-tests via `anova()`) are extracted:

- **Global LRT** (reduced vs. full) — tests whether TSS explains exon inclusion beyond time alone.
- **Interaction LRT** (full vs. interaction) — tests whether the effect of time differs across TSSs (i.e., whether the relationship between exon inclusion and TSS usage changes across time).

Estimated marginal means (EMMs) are computed from the interaction model using `emmeans(~ Time | TSS)`, providing per-site, per-timepoint PSI estimates on the response scale with confidence intervals.

Convergence is monitored: any event where the interaction model fails to converge (either via the `$converged` flag or a detected `"algorithm did not converge"` warning) is flagged and excluded from the final outputs.

---

## 6. FDR correction

P-values from both LRTs are extracted at event level (one value per `event_id`), and Benjamini–Hochberg FDR correction is applied across all events separately for each test. Corrected q-values are joined back to the full coefficient table.

---

## Outputs

| File | Description |
|---|---|
| `all_events_quasi_TSS.csv` | One row per model term per event. Contains raw and FDR-corrected p-values for both LRTs, convergence flag, and model family. |
| `all_events_emmeans_quasi_TSS.csv` | One row per `(event_id, TSS_id, Time)` combination. Contains EMM-derived PSI estimates and confidence intervals. Used for visualisation and effect-size calculation. |

---
