from src.ryp import r, to_r
import polars as pl

r(
    """
    library(dplyr)
    library(ggplot2)

    colorVector <- c(
        "full-splice_match" = "#009E73",
        "incomplete-splice_match"   = "#0072B2",
        "novel_in_catalog"   = "#D55E00",
        "novel_not_in_catalog"   = "#E69F00",
        "Other" = "#000000"
    )

    colorVector <- c(
        "FSM" = "#009E73",
        "ISM"   = "#0072B2",
        "NIC"   = "#D55E00",
        "NNC"   = "#E69F00",
        "Other" = "#000000"
    )
    """
)

def plot_transcript_class_hist(classification):
    structural_category_labels = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM", "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC", "Other": "Other"}
    classification = classification\
        .with_columns(
            structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
                .then(pl.col("structural_category"))\
                .otherwise(pl.lit("Other"))
        )\
        .with_columns(
            structural_category2 = pl.col("structural_category2").map_elements(lambda x: structural_category_labels[x])
        )

    summary_df = classification\
        .group_by("structural_category2")\
        .agg(pl.len().alias("len"))\
        .with_columns(
            (pl.col("len") / pl.col("len").sum() * 100).alias("percentage")
        )
    
    to_r(summary_df, "summary_df")

    r(
        """
        summary_df$structural_category2 <- factor(summary_df$structural_category2, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))

        summary_df %>% 
            ggplot(aes(x = structural_category2, y = len, fill = structural_category2)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "white", size = 5) +
            theme(axis.text.x = element_text(size = 18)) +
            labs(x = "Structural Category", y = "Percentage", title = "Transcripts Identified by Novelty") +
            scale_fill_manual("Structural Category", values=colorVector) +
            scale_y_continuous(labels = function(x) x/1000) +
            ylab("Transcripts (x 10^3)") +
            xlab(NULL)

        ggsave("transcript_classification_hist.png", width=7, height=4)
        """
    )