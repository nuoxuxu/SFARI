# This script takes gff3 made by TransDecoder and shapes it into a format
# that can be used to construct Annotation by GenomeKit

from src.utils import read_gff, read_gtf
import polars as pl
from src.single_cell import SingleCell

# Read TD_gff and clean

TD_gff_path = "full_nt.fasta.transdecoder.genome_remove_spaces.gff3"

TD_gff = read_gff(TD_gff_path, attributes=["ID", "Parent"])\
    .with_columns(
        pl.col("feature").replace({"gene": "gene", "five_prime_UTR": "five_prime_UTR", "exon": "exon", "mRNA": "transcript", "CDS": "CDS", "three_prime_UTR": "three_prime_UTR"})
    ).drop_nulls("feature")

# Format ID and Parent columns and add exon_number column

TD_gff = TD_gff\
    .with_columns(
        pl.when(pl.col("feature")=="gene")\
            .then(pl.col("ID").str.split("^").map_elements(lambda x: x[0], return_dtype=pl.String))\
            .when(pl.col("feature")=="transcript")\
            .then(pl.col("ID").str.extract(r"^(.*)\.p\d+$"))\
            .otherwise(pl.col("ID")),
        pl.when(pl.col("feature")=="transcript")\
            .then(pl.col("Parent").str.extract(r"^(.*)\^chr"))\
            .when(pl.col("feature").is_in(["three_prime_UTR", "five_prime_UTR", "CDS", "exon"]))\
            .then(pl.col("Parent").str.extract(r"^(.*)\.p\d+$"))\
            .otherwise(pl.col("Parent")),
        exon_number = pl.when(pl.col("feature")=="exon")\
            .then(pl.col("ID").str.extract(r"exon(\d+)"))\
            .otherwise(0)
    )

# Add gene_id and transcript_id

TD_gff = TD_gff\
    .with_columns(
        gene_id = pl.when(pl.col("feature")=="gene")\
        .then(pl.col("ID"))\
        .when(pl.col("feature")=="transcript")\
        .then(pl.col("Parent"))\
        .otherwise(pl.col("Parent").str.split(".").map_elements(lambda x: ".".join([x[0], x[1]]), return_dtype=pl.String)),

        transcript_id = pl.when(pl.col("feature")=="gene")\
        .then(pl.lit(None))\
        .when(pl.col("feature")=="transcript")\
        .then(pl.col("ID"))\
        .otherwise(pl.col("Parent"))
    )

# TD_gff\
#     .with_columns(
#         attributes = pl.when(pl.col("feature")=="gene")\
#             .then(pl.lit("ID=")+pl.col("ID"))\
#             .when(pl.col("feature").is_in(["transcript", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
#             .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
#             .when(pl.col("feature")=="exon")\
#             .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number"))
#     ).drop(["ID", "Parent", "exon_number"])\
#     .write_csv("full_nt.fasta.transdecoder.genome.genomekit.gff3", include_header=False, quote_style="never", separator="\t")

path_to_merged_collapsed_gff = "proc/merged_collapsed.filtered.gff"
merged_collapsed_gff = read_gtf(path_to_merged_collapsed_gff, keep_attributes=False, attributes=["gene_id", "transcript_id"])
gene_gtf = merged_collapsed_gff\
    .filter(pl.col("feature") == "transcript")\
    .group_by("gene_id", maintain_order=True)\
    .agg(
        pl.col("seqname").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        pl.col("source").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        pl.col("feature").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        pl.col("start").min(),
        pl.col("end").max(),
        pl.col("score").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        pl.col("strand").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        pl.col("frame").unique().map_elements(lambda x: x[0], return_dtype=pl.String)
    )\
    .with_columns(
        pl.lit(None).alias("transcript_id"),
        pl.col("feature").replace({"transcript": "gene"})
    )\
    .with_columns(
        transcript_id = pl.col("transcript_id").cast(pl.String)
    )['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id']

PB_gtf = pl.concat([gene_gtf, merged_collapsed_gff], how="vertical")
#     .with_columns(
#         pl.col("feature").cast(pl.Enum(["gene", "transcript", "exon"]))
#     )

# combined_gtf = combined_gtf.sort("gene_id", "transcript_id", "feature", "start")

# Add exon_number

exon_number = PB_gtf\
    .filter(pl.col("feature")=="exon")\
    .group_by("transcript_id")\
    .agg(
        pl.col("start"),
        pl.cum_count("start").alias("exon_number")
    ).explode("start", "exon_number")
PB_gtf = PB_gtf\
    .join(exon_number, on=["transcript_id", "start"], how="left")

# Add ID and Parent

PB_gtf = PB_gtf\
    .with_columns(
        ID = pl.when(pl.col("feature")=="gene")\
            .then(pl.col("gene_id"))\
            .when(pl.col("feature")=="transcript")\
            .then(pl.col("transcript_id"))\
            .when(pl.col("feature")=="exon")\
            .then(pl.col("transcript_id")+pl.lit(".exon")+pl.col("exon_number").cast(pl.String)),
        Parent = pl.when(pl.col("feature")=="gene")\
            .then(None)\
            .when(pl.col("feature")=="transcript")\
            .then(pl.col("gene_id"))\
            .when(pl.col("feature")=="exon")\
            .then(pl.col("transcript_id"))
    )

# Keep only pbids that are not in TD_gff and their associated gene entries

lr_bulk = SingleCell("results/long_read/pbid_filtered.h5ad")
pbid_not_in_TD_gff = set(lr_bulk.var["pbid"].to_list()) - set(TD_gff.filter(pl.col("feature")=="transcript").unique("ID")["ID"].to_list())
genes_for_pbid_not_in_TD_gff = list(set([(lambda x: ".".join([x.split(".")[0], x.split(".")[1]]))(pbid) for pbid in pbid_not_in_TD_gff]))
PB_gtf = PB_gtf\
    .filter((pl.col("feature")=="gene").and_(pl.col("ID").is_in(genes_for_pbid_not_in_TD_gff))|((pl.col("feature")!="gene").and_(pl.col("transcript_id").is_in(pbid_not_in_TD_gff))))

# sort PB_gtf

PB_gtf = PB_gtf\
    .with_columns(
        feature = pl.col("feature").cast(pl.Enum(["gene", "transcript", "exon"]))
    )\
    .sort("gene_id", "transcript_id", "feature", "start")

# Combine PB_gtf with TD_gff

combined_gtf = pl.concat(
    [PB_gtf['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number', 'ID', 'Parent']\
        .with_columns(
            exon_number=pl.col("exon_number").cast(pl.String),
            feature = pl.col("feature").cast(pl.String)
        ), 
    TD_gff['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'exon_number', 'ID', 'Parent']],
    how="vertical"
)

# Remove duplicated gene entries

duplicate_genes = combined_gtf.filter(pl.col("feature")=="gene").group_by("ID").count().filter(pl.col("count")==2)
combined_gtf = combined_gtf\
    .filter(((pl.col("feature")=="gene").and_(pl.col("source")=="transdecoder").and_(pl.col("ID").is_in(duplicate_genes["ID"]))).not_())

# Sort combined_gtf so that features of the same gene_id are grouped together followed then by transcript_id

combined_gtf = combined_gtf\
    .sort("gene_id", "transcript_id", maintain_order=True)

# Get Attributes and drop redundant columns

combined_gtf = combined_gtf\
    .with_columns(
        attributes = pl.when(pl.col("feature")=="gene")\
            .then(pl.lit("ID=")+pl.col("ID"))\
            .when(pl.col("feature").is_in(["transcript", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
            .when(pl.col("feature")=="exon")\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number").cast(pl.String))
    )\
    .drop("exon_number", "ID", "Parent", "gene_id", "transcript_id")\
    ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

combined_gtf.write_csv("complete_TD_genomekit.gff3", include_header=False, quote_style="never", separator="\t")

combined_gtf\
    .filter(
        pl.col("source") == "PacBio"
    ).write_csv("test.gff3", include_header=False, quote_style="never", separator="\t")
#----------------------------------------------------------------------
PB_gtf = PB_gtf\
    .filter(pl.col("transcript_id").is_in(pbid_not_in_TD_gff).or_((pl.col("feature")=="gene").and_(pl.col("gene_id").is_in(gene_id_not_in_TD_gff))))\
    .with_columns(
        attributes = pl.when(pl.col("feature")=="gene")\
            .then(pl.lit("ID=")+pl.col("ID"))\
            .when(pl.col("feature").is_in(["transcript", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
            .when(pl.col("feature")=="exon")\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number").cast(pl.String))
    ).drop("exon_number", "ID", "Parent", "gene_id", "transcript_id")

TD_gff = TD_gff\
    .with_columns(
        attributes = pl.when(pl.col("feature")=="gene")\
            .then(pl.lit("ID=")+pl.col("ID"))\
            .when(pl.col("feature").is_in(["transcript", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
            .when(pl.col("feature")=="exon")\
            .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number"))
    ).drop(["ID", "Parent", "exon_number"])\
    .with_columns(
        feature = pl.col("feature").replace({"gene": "gene", "transcript": "transcript", "five_prime_UTR": "five_prime_UTR", "three_prime_UTR": "three_prime_UTR", "exon": "exon", "CDS": "exon"}),
        feature_old = pl.col("feature")
    )

pl.concat([PB_gtf.with_columns(pl.col("feature").cast(pl.String)), TD_gff], how="vertical")\
    .with_columns(
        pl.col("feature").cast(pl.Enum(["gene", "transcript", "five_prime_UTR", "exon", "three_prime_UTR"]))
    )\
    .sort("feature")