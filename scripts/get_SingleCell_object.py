import polars as pl
from src.single_cell import SingleCell
from scipy.sparse import csr_array
import argparse

def main():
    parser = argparse.ArgumentParser(description='Get single_cell object from long read data')
    parser.add_argument('--id_to_sample', action='store', dest='id_to_sample', type=str, required=True)
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--read_stat', action='store', dest='read_stat', type=str, required=True)
    parser.add_argument('--output', action='store', dest='output', type=str, required=True)
    params = parser.parse_args()

    id_to_sample = pl.read_csv(params.id_to_sample, separator = "\t", has_header = False, new_columns = ["id", "sample"])\
        .with_columns(
            pl.col("id").map_elements(lambda s: s.split("/")[0], return_dtype = pl.String),
            pl.col("sample").map_elements(lambda s: s.rsplit("/")[3].rsplit("_", 2)[0], return_dtype = pl.String)
        )\
        .to_pandas().set_index("id").to_dict()["sample"]

    categories_to_show = ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]
    
    classification = pl.read_csv(params.classification, separator="\t", null_values=["NA"])\
        .rename({"isoform": "pbid"})\
        .with_columns(
            pl.when(pl.col("structural_category").is_in(categories_to_show)).then(pl.col("structural_category")).otherwise(pl.lit("Other")).alias("structural_category2")
            )

    read_stat = pl.read_csv(params.read_stat, separator = "\t")\
        .with_columns(
            pl.col("id").map_elements(lambda s: s.split("/")[0], return_dtype = pl.String),pl.lit(1).alias("count")
        )\
        .group_by("id", "pbid")\
        .agg(pl.sum("count"))\
        .sort("pbid", descending=False)\
        .pivot("id", index = "pbid", values = "count")\
        .fill_null(0)\
        .rename(id_to_sample)\
        ["pbid", "iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2"]

    sampleID = read_stat.columns[1:]
    time_point = [x.split("_")[0] for x in sampleID]
    myDesign = pl.DataFrame(
        {
            "sampleID": sampleID,
            "time_point": time_point
        }
    )

    # Matching pbid in classification and read_stat

    read_stat = classification\
        .join(read_stat, on = "pbid", how = "left")\
        .select(read_stat.columns)

    classification = classification\
        .join(read_stat, on = "pbid", how = "left")\
        .select(classification.columns)

    # Create SingleCell object

    lr_bulk = SingleCell(
        X = csr_array(read_stat.drop("pbid").to_numpy().T),
        obs = myDesign,
        var = classification
    )

    # Convert columns that contain missing values only to float64

    lr_bulk.var = lr_bulk.var\
        .with_columns(
            pl.col(col).cast(pl.Float64) for col in lr_bulk.var.columns if (lr_bulk.var.select(pl.col(col).null_count() == 4897888)).to_numpy()[0]
        )

    # Save SingleCell object
    lr_bulk.save(params.output)