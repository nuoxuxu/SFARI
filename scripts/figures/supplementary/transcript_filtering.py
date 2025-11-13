from ryp import r, to_r, to_py
import polars as pl
import polars.selectors as cs

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
            )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
            ).drop("attributes")

def get_classification(classification_path):
    print("Reading isoform classification file...")
    classification = pl.read_csv(classification_path, separator="\t", null_values=["NA"])
    print(f"{classification.shape[0]} isoforms are in classification file, these are the isoforms that passes pigeon filter")
    return classification

def get_exp_tx(expression, min_reads=5, min_n_sample=2):
    return expression\
        .with_columns(
            cs.numeric() > min_reads
        )\
        .filter(
            pl.sum_horizontal(cs.boolean()) > min_n_sample
        )\
        ["isoform"].to_list()

def get_polyA_tx(gtf, polyA_site):
    validated_pbids = gtf\
        .filter(pl.col("feature")=="transcript")\
        .select(
            pl.col("seqname"),
            pl.col("transcript_id"),
            pos = pl.when(pl.col("strand")=="+")\
                .then(pl.col("end"))\
                .otherwise(pl.col("start"))
        )\
        .join_where(
            polyA_site,
            (pl.col("pos") >= pl.col("start")) &
            (pl.col("pos") <= pl.col("end"))
        )\
        .filter(
            pl.col("seqname")==pl.col("chrom")
        )\
        ["transcript_id"].to_list()
    return validated_pbids
    
def get_CAGE_tx(gtf, reftss):
    validated_pbids = gtf\
        .filter(pl.col("feature")=="transcript")\
        .select(
            pl.col("seqname"),
            pl.col("transcript_id"),
            pos = pl.when(pl.col("strand")=="+")\
                .then(pl.col("start"))\
                .otherwise(pl.col("end"))
        )\
        .join_where(
            reftss,
            (pl.col("pos") >= (pl.col("start")-100)) &
            (pl.col("pos") <= (pl.col("end")+100))
        )\
        .filter(
            (pl.col("seqname") == pl.col("chrom"))
        )\
        .unique("transcript_id")\
        ["transcript_id"].to_list()
    return validated_pbids

# Get pigeon filtered classification
classification = get_classification("nextflow_results/V47/merged_collapsed_classification.filtered_lite_classification.txt")
filtered_gff = read_gtf("nextflow_results/V47/merged_collapsed.sorted.filtered_lite.gff")

# Remove lowly expressed isoforms
expression = pl.read_parquet("nextflow_results/V47/full_expression.parquet")
expression = expression\
    .filter(
        pl.col("isoform").is_in(classification["isoform"])
    )
expression = expression\
    .filter(
        pl.col("isoform").is_in(get_exp_tx(expression, 5, 2))
    )
classification = classification\
    .filter(
        pl.col("isoform").is_in(expression["isoform"])
    )

# Add CAGE and polyA site information to classification
reftss = pl.read_csv(
        "./data/refTSS_v3.3_human_coordinate.hg38.sorted.bed", 
        separator="\t", has_header = False, 
        new_columns=["chrom", "start", "end", "name", "score", "strand"]
    )
polyA_site = pl.read_csv(
        "./data/atlas.clusters.2.0.GRCh38.96.bed", 
        separator="\t", new_columns=["chrom", "start", "end", "name", "score", "strand"], 
        schema_overrides={"chrom": pl.String}
    )

classification = classification\
    .with_columns(
        within_CAGE_peak = pl.col("isoform").is_in(get_CAGE_tx(filtered_gff, reftss)),
        within_polyA_site = pl.col("isoform").is_in(get_polyA_tx(filtered_gff, polyA_site))
    )

# Add new column structural_category2
classification = classification\
    .with_columns(
        structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
            .then(pl.col("structural_category"))\
            .otherwise(pl.lit("Other"))
    )\
    .with_columns(
        structural_category2 = pl.when(pl.col("subcategory") == "3prime_fragment")\
            .then(pl.lit("3prime_fragment"))\
            .when(pl.col("subcategory") == "5prime_fragment")\
            .then(pl.lit("5prime_fragment"))\
            .when(pl.col("subcategory") == "internal_fragment")\
            .then(pl.lit("internal_fragment"))\
            .when(pl.col("subcategory") == "intron_retention")\
            .then(pl.lit("intron_retention"))\
            .otherwise(pl.col("structural_category2"))
    )

classification = classification\
    .with_columns(
        pl.col("structural_category2").replace({"incomplete-splice_match": "ISM", 
                                               "full-splice_match": "FSM", 
                                               "novel_in_catalog": "NIC", 
                                               "novel_not_in_catalog": "NNC", 
                                               "Other": "Other"}).alias("structural_category2")
    )

# Calculate polyA site percentage per structural category and transfer to R object
structural_category_polyA_site_percentage = classification\
    .group_by(["structural_category2", "within_polyA_site"])\
    .len()\
    .group_by("structural_category2")\
    .agg([
        pl.col("len").filter(pl.col("within_polyA_site") == True).sum().alias("true_len"),
        pl.col("len").sum().alias("total_len")
    ])\
    .with_columns(
        (pl.col("true_len") / pl.col("total_len") * 100).alias("polyA_site_percentage")
    )
to_r(structural_category_polyA_site_percentage, "structural_category_polyA_site_percentage")

# Calculate CAGE peak percentage per structural category and transfer to R object
structural_category_CAGE_peak_percentage = classification\
    .group_by(["structural_category2", "within_CAGE_peak"])\
    .len()\
    .group_by("structural_category2")\
    .agg([
        pl.col("len").filter(pl.col("within_CAGE_peak") == True).sum().alias("true_len"),
        pl.col("len").sum().alias("total_len")
    ])\
    .with_columns(
        (pl.col("true_len") / pl.col("total_len") * 100).alias("CAGE_percentage")
    )
to_r(structural_category_CAGE_peak_percentage, "structural_category_CAGE_peak_percentage")

r(
    """
    library(dplyr)
    library(ggplot2)
    library(arrow)
    library(patchwork)

    colorVector <- c(
        "FSM" = "#009E73",
        "ISM" = "#0072B2",
        "NIC" = "#D55E00",
        "NNC" = "#E69F00",
        "Other" = "#000000"
    )

    structural_category_labels <- c(
        "full-splice_match"        = "FSM",
        "incomplete-splice_match"  = "ISM",
        "novel_in_catalog"         = "NIC",
        "novel_not_in_catalog"     = "NNC",
        "Other"                    = "Other"
    )
    my_theme <- theme_classic() +
        theme(
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            legend.position = "none"
        )
    update_geom_defaults("text", list(size = 5))    
    theme_set(my_theme)
    """
)

r(
    """
    p1 <- structural_category_polyA_site_percentage %>% 
        mutate(
            structural_category2 = factor(structural_category2, levels = c("FSM", "5prime_fragment", "3prime_fragment", "intron_retention", "internal_fragment", "NIC", "NNC", "Other"))
        ) %>%
        filter(structural_category2 != "incomplete-splice_match") %>% 
        ggplot(aes(x = structural_category2, y = polyA_site_percentage)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(polyA_site_percentage, 1), "%")), hjust = 2, colour = "black", size = 3.5) +
        theme_minimal() +
        labs(x = "Structural Category", y = "% within polyA site") +
        coord_flip()

    p2 <- structural_category_CAGE_peak_percentage %>% 
        mutate(
            structural_category2 = factor(structural_category2, levels = c("FSM", "5prime_fragment", "3prime_fragment", "intron_retention", "internal_fragment", "NIC", "NNC", "Other"))
        ) %>%    
        filter(structural_category2 != "incomplete-splice_match") %>% 
        ggplot(aes(x = structural_category2, y = CAGE_percentage)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(CAGE_percentage, 1), "%")), hjust = 2, colour = "black", size = 3.5) +
        theme_minimal() +
        labs(x = "Structural Category", y = "% within CAGE peak") +
        coord_flip()        
    p1 + p2 + plot_layout(ncol=2)
    ggsave("figures/supplementary/ISM_filtering.pdf", width=7.5, height=2.5)
    """
)
