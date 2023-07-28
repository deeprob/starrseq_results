import os
import numpy as np
import pandas as pd
from itertools import combinations


def create_homer_bed(df, libn1, libn2, save_file):
    if len(df)==0:
        return
    # drop duplicates
    df = df.fillna("NA").drop_duplicates(keep="first")
    # split the merged coordinates to bed format
    df = pd.concat((df, df.chrom_coord.str.split("_", expand=True).rename(columns={0: "chrom", 1: "start", 2: "end"})), axis=1)
    # add strand info
    df["strand"] = "."
    df = df.loc[:, ["chrom", "start", "end", "chrom_coord", "CC_peak", "strand", f"{libn1}_log2FoldChange_act", f"{libn2}_log2FoldChange_act", f"gene_name", f"{libn1}_log2FoldChange_exp", f"{libn2}_log2FoldChange_exp", f"{libn1}_peak", f"{libn2}_peak"]]
    return df.to_csv(save_file, sep="\t", header=False, index=False)


def read_and_filter_df(filename, libraries):
    df = pd.read_csv(filename)
    df = df.dropna()
    # proportional query
    query = " & ".join([
        f"((`{lib}_log2FoldChange_act`<0) & (`{lib}_log2FoldChange_exp`<0) | (`{lib}_log2FoldChange_act`>0) & (`{lib}_log2FoldChange_exp`>0))" for lib in libraries])
    df = df.query(query)
    return df


if __name__=="__main__":
    abc_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks/abc_da_de_table_peaks.csv"
    nearest_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks/nearest_da_de_table_peaks.csv"
    lib_names = ["ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2"]
    save_dir = os.path.join("/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks", "discordant")

    combos = list(combinations(lib_names,2))
    discordant_counts = []
    for lib_combo in combos:
        # combine abc and nearest df
        abc_df = read_and_filter_df(abc_file, lib_combo)
        nearest_df = read_and_filter_df(nearest_file, lib_combo)
        df = pd.concat((abc_df, nearest_df)).drop_duplicates(subset=["chrom_coord", "gene_name"])
        df = df.loc[((df["CC_peak"]==1|(df[f"{lib_combo[0]}_peak"]==1)|(df[f"{lib_combo[1]}_peak"]==1)))]
        df["eg_pairs"] = df.chrom_coord + "::" + df.gene_name
        df = df.set_index("eg_pairs")
        # query for induced in one repressed in another or vice versa
        query = f"(`{lib_combo[0]}_log2FoldChange_act`>0) & (`{lib_combo[1]}_log2FoldChange_act`<0) | (`{lib_combo[0]}_log2FoldChange_act`<0) & (`{lib_combo[1]}_log2FoldChange_act`>0)"
        df = df.query(query)
        df = df.loc[:, ["chrom_coord", "CC_peak", "gene_name"] + [c for c in df.columns if c.startswith(lib_combo)]]
        save_file = os.path.join(save_dir, f"{lib_combo[0]}_vs_{lib_combo[1]}", "discordant.bed")
        os.makedirs(os.path.dirname(save_file), exist_ok=True)
        create_homer_bed(df, lib_combo[0], lib_combo[1], save_file)
