from src.ryp import r, to_r
import polars as pl

r(
    """
    library(dplyr)
    library(ggplot2)
    library(scales)

    colorVector <- c(
        "FSM" = "#009E73",
        "ISM"   = "#0072B2",
        "NIC"   = "#D55E00",
        "NNC"   = "#E69F00",
        "Other" = "#000000"
    )

    proteinColorVector <- c(
        "pFSM" = "#009E73",
        "pISM"   = "#0072B2",
        "pNIC"   = "#D55E00",
        "pNNC"   = "#E69F00",
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
            theme_bw() +
            theme(
                axis.text.x = element_text(size = 18),
                panel.grid.minor = element_blank()
                ) +
            labs(x = "Structural Category", y = "Percentage", title = "Transcripts Identified by Novelty") +
            scale_fill_manual("Structural Category", values=colorVector) +
            scale_y_continuous(labels = function(x) x/1000) +
            ylab("Transcripts (x 10^3)") +
            xlab(NULL)

        ggsave("transcript_classification_hist.png", width=7, height=4)
        """
    )

def plot_protein_class_hist(protein_classification, orfanage_gtf, save_path):
    from src.utils import collapse_isoforms_to_proteoforms
    
    isoforms_to_proteoforms = collapse_isoforms_to_proteoforms(orfanage_gtf)

    protein_classification = protein_classification\
        .rename({"pb": "isoform"})\
        .join(isoforms_to_proteoforms, on="isoform")\
        .unique("base_isoform")\
        .with_columns(
            structural_category2 = pl.when(pl.col("protein_classification_base").is_in(["pNIC", "pFSM", "pISM", "pNNC"]))\
                .then(pl.col("protein_classification_base"))\
                .otherwise(pl.lit("Other"))
        )

    summary_df = protein_classification\
        .group_by("structural_category2")\
        .agg(pl.len().alias("len"))\
        .with_columns(
            (pl.col("len") / pl.col("len").sum() * 100).alias("percentage")
        )

    to_r(summary_df, "summary_df")
    to_r(save_path, "save_path")

    r(
        """
        library(ggplot2)
        library(dplyr)
        colorVector <- c(
            "pFSM" = "#009E73",
            "pISM"   = "#0072B2",
            "pNIC"   = "#D55E00",
            "pNNC"   = "#E69F00",
            "Other" = "#000000"
        )

        summary_df$structural_category2 <- factor(summary_df$structural_category2, levels = c("pFSM", "pISM", "pNIC", "pNNC", "Other"))
        summary_df %>%
            filter(structural_category2 %in% c("pFSM", "pNIC", "pNNC")) %>% 
            ggplot(aes(x = structural_category2, y = len, fill = structural_category2)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "white", size = 5) +
            theme_bw() +
            theme(
                axis.text.x = element_text(size = 18),
                panel.grid.minor = element_blank()
                ) +
            labs(x = "Structural Category", y = "Percentage", title = "Proteoforms Identified by Novelty") +
            scale_fill_manual("Structural Category", values=colorVector) +
            scale_y_continuous(labels = function(x) x/1000) +
            ylab("Proteoforms (x 10^3)") +
            xlab(NULL)

        ggsave(save_path, width=7, height=6)
        """
    )

def plot_nexon_vs_abundance(classification, expression):
    structural_category_labels = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM", "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC", "Other": "Other"}
    classification = classification\
        .with_columns(
            structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
                .then(pl.col("structural_category"))\
                .otherwise(pl.lit("Other"))
        )\
        .with_columns(
            structural_category2 = pl.col("structural_category2").replace_strict(structural_category_labels)
        )
    df = pl.concat(
        [classification["isoform", "structural_category2", "exons"], expression[:, 1:].sum_horizontal().to_frame()],
        how = "horizontal"
        )
    to_r(df, "df")
    r(
        """
            df %>%
                ggplot(aes(x=exons, fill=structural_category2)) +
                geom_histogram(alpha=.75, binwidth = 1) +
                theme_bw() + 
                xlim(1,40) + 
                scale_fill_manual(values=colorVector) + 
                labs(x="# Exons", y="Transcripts(x 10^3)") + 
                scale_fill_manual("Structural Category", values=colorVector) +
                ggtitle('Exons per  Transcript') + 
                scale_y_continuous(labels = function(x) x/1000)
        """
    )

def plot_len_vs_abundance(classification, expression):
    structural_category_labels = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM", "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC", "Other": "Other"}
    classification = classification\
        .with_columns(
            structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
                .then(pl.col("structural_category"))\
                .otherwise(pl.lit("Other"))
        )\
        .with_columns(
            structural_category2 = pl.col("structural_category2").replace_strict(structural_category_labels)
        )
    df = pl.concat(
        [classification["isoform", "structural_category2", "length"], expression[:, 1:].sum_horizontal().to_frame()],
        how = "horizontal"
        )
    to_r(df, "df")
    r(
        """
            df %>%
                ggplot(
                    aes(x=length, fill=structural_category2)
                    ) +
                geom_histogram(alpha=.75) +
                theme_bw() + 
                scale_fill_manual(values=colorVector) +
                scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)),limits = c(300,10^4)) +
                scale_y_continuous(labels = function(x) x/1000) +
                ylab("Transcripts (x 10^3)") +
                labs(x="Transcript Length (bp)") +
                scale_fill_manual("Structural Category", values=colorVector) +
                ggtitle("Transcript length distribution")
        """
    )

  

def plot_transcript_class_abundance(classification, expression):
    structural_category_labels = {"full-splice_match": "FSM", "incomplete-splice_match": "ISM", "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC", "Other": "Other"}
    classification = classification\
        .with_columns(
            structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
                .then(pl.col("structural_category"))\
                .otherwise(pl.lit("Other"))
        )\
        .with_columns(
            structural_category2 = pl.col("structural_category2").replace_strict(structural_category_labels)
        )
    df = pl.concat(
        [classification["isoform", "structural_category2"], expression[:, 1:].sum_horizontal().to_frame()],
        how = "horizontal"
        )
    to_r(df, "df")
    r(
        """
        library(ggplot2)
        library(dplyr)
        df %>% 
            filter(sum > 10) %>%
            ggplot(aes(x=sum, fill=structural_category2)) +
            geom_histogram(position=position_fill(), alpha=.75) +
            theme_bw() +
            scale_x_log10() +
            annotation_logticks(scaled = T, sides= "b") +
            theme(panel.grid.minor = element_blank()) + 
            labs(x="Min observed counts", y="Transcript proportion") + 
            ggtitle("Transcript novelty vs Abundance") + 
            theme(plot.title = element_text(hjust=.5)) + 
            scale_fill_manual("Structural Category", values=colorVector)    

        """
    )