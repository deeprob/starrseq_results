import os
import numpy as np
from scipy.stats import pearsonr
import pandas as pd


def get_active_enhancer_links(activity_file, libraries, save_file):
    activity_df = pd.read_csv(activity_file, low_memory=False)
    # select only active enhancers in any library
    activity_query = " | ".join([f"(`{lib}_peak`==1)" for lib in libraries])
    activity_df = activity_df.query(activity_query)
    # remove low expression
    expression_query = " | ".join([f"(`{lib}_exp`>=1)" for lib in libraries])
    activity_df = activity_df.query(expression_query)
    # # keep their target gene column, target gene exp columns and activity columns and chrom_coord column
    cols_to_keep = ["chrom_coord", "gene_name", "per_gene_corr"]
    cols_to_keep = cols_to_keep + [f"{lib}_act" for lib in libraries]
    cols_to_keep = cols_to_keep + [f"{lib}_peak" for lib in libraries]
    cols_to_keep = cols_to_keep + [f"{lib}_log2FoldChange_act" for lib in libraries[1:]]
    cols_to_keep = cols_to_keep + [f"{lib}_padj_act" for lib in libraries[1:]]
    cols_to_keep = cols_to_keep + [f"{lib}_exp" for lib in libraries]
    cols_to_keep = cols_to_keep + [f"{lib}_log2FoldChange_exp" for lib in libraries[1:]]
    cols_to_keep = cols_to_keep + [f"{lib}_padj_exp" for lib in libraries[1:]]
    activity_df = activity_df.loc[~activity_df.per_gene_corr.isna(), cols_to_keep]
    activity_df.to_csv(save_file, index=False)
    return activity_df

def get_active_enhancer_lib_links(activity_df, lib_name, cc_name):
    # select only active enhancers
    activity_query = " | ".join([f"(`{lib}_peak`==1)" for lib in [cc_name, lib_name]])
    activity_df = activity_df.query(activity_query)
    # remove low expression
    expression_query = " | ".join([f"(`{lib}_exp`>=1)" for lib in [cc_name, lib_name]])
    activity_df = activity_df.query(expression_query)
    # select sda and sde only
    da_de_query = f"(`{lib_name}_padj_act`<=0.05) & (`{lib_name}_padj_exp`<=0.05)"
    activity_df = activity_df.query(da_de_query)
    # filter for specific columns
    cols_to_keep = [
        "chrom_coord", "gene_name", f"{cc_name}_act", f"{cc_name}_exp", f"{cc_name}_peak",
        f"{lib_name}_log2FoldChange_act", f"{lib_name}_log2FoldChange_exp", f"{lib_name}_peak",
        f"{lib_name}_act", f"{lib_name}_padj_act", f"{lib_name}_padj_exp"
        ]
    if "distance" in activity_df.columns:
        cols_to_keep += ["distance"]
    activity_df = activity_df.loc[:, cols_to_keep]
    return activity_df

def merge_multi_fragment_into_one(merge_df, lib):
    # merge overlapping locations
    # sort the dataframe by chromosomal coordinates
    merge_df_locations = merge_df.chrom_coord.str.split("_", expand=True)
    merge_df_locations.columns = ["chr", "start", "end"]
    merge_df_locations = merge_df_locations.astype({'start': 'int32', 'end': 'int32'})
    merge_df = merge_df.merge(merge_df_locations, left_index=True, right_index=True)
    merge_df = merge_df.sort_values(["chr", "start", "end"])

    dict_function = {
        "chrom_coord": lambda x: list(x)[0],
        "gene_name": lambda x: list(x), 
        "CC_act": max, # max CC activity
        f"{lib}_act": max, # max activity value for each fragment
        f"{lib}_padj_act": min, # min significance value for each fragment
        f"{lib}_log2FoldChange_act": lambda x: max(x) if all(i>0 for i in x) else min(x),
        f"{lib}_padj_exp": min, # min significance value for each fragment
        f"{lib}_log2FoldChange_exp": lambda x: max(x) if all(i>0 for i in x) else min(x),
        "CC_peak": max, 
        f"{lib}_peak": max,
        "distance": min
        }
        
    grouped_dfs = []
    for group, df in merge_df.groupby(["chr", "gene_name"]): # https://stackoverflow.com/questions/46732760/merge-rows-pandas-dataframe-based-on-condition
        grouped_df = df.groupby(((df.start  - df.end.shift(1)) > 0).cumsum()).agg({k:v for k,v in dict_function.items() if k in merge_df.columns})
        grouped_df["chr"] = group[0]
        grouped_df["gene_name"] = group[1]
        grouped_dfs.append(grouped_df)

    merged_df = pd.concat(grouped_dfs)
    return merged_df.drop(columns="chr")


if __name__ == "__main__":
    abc_activity_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/abc_da_de_table.csv"
    nearest_activity_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/nearest_da_de_table.csv"
    nearest_target_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/predictions/nearest_genes/closest.bed"

    abc_save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks/abc_da_de_table_peaks.csv"
    nearest_save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks/nearest_da_de_table_peaks.csv"
    
    lib_links_save_dir = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks"
    libraries = ["CC", "ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]

    # For each method get the da de table for active enhancers only
    abc_df = get_active_enhancer_links(abc_activity_file, libraries, abc_save_file)
    nearest_df = get_active_enhancer_links(nearest_activity_file, libraries, nearest_save_file)
    nearest_target_df = pd.read_csv(
        nearest_target_file, sep="\t", header=None, 
        names=["chrm", "start", "end", "gene_name", "strand", "distance"], 
        usecols=[0,1,2,7,8,9]
        )
    nearest_target_df["chrom_coord"] = nearest_target_df.chrm + "_" + nearest_target_df.start.astype(str) + "_" + nearest_target_df.end.astype(str)
    nearest_df = nearest_df.merge(nearest_target_df, on=["chrom_coord", "gene_name"])


    for lib_name in libraries[1:]:
        # get sda sde enhancers per library
        lib_df = get_active_enhancer_lib_links(abc_df, lib_name, "CC")
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "abc_sda_sde_table_peaks_fragments.csv")
        os.makedirs(os.path.dirname(lib_save_file), exist_ok=True)
        lib_df.to_csv(lib_save_file, index=False)
        lib_df = merge_multi_fragment_into_one(lib_df, lib_name)
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "abc_sda_sde_table_peaks.csv")
        os.makedirs(os.path.dirname(lib_save_file), exist_ok=True)
        lib_df.to_csv(lib_save_file, index=False)
        # find proportional enhancers-gene pairs per library
        lib_df = lib_df.loc[np.sign(lib_df[f"{lib_name}_log2FoldChange_act"])==np.sign(lib_df[f"{lib_name}_log2FoldChange_exp"])]
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "abc_sda_sde_table_peaks_proportional.csv")
        os.makedirs(os.path.dirname(lib_save_file), exist_ok=True)
        lib_df.to_csv(lib_save_file, index=False)

        lib_df = get_active_enhancer_lib_links(nearest_df, lib_name, "CC")
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "nearest_sda_sde_table_peaks_fragments.csv")
        lib_df.to_csv(lib_save_file, index=False)        
        lib_df = merge_multi_fragment_into_one(lib_df, lib_name)
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "nearest_sda_sde_table_peaks.csv")
        lib_df.to_csv(lib_save_file, index=False)
        # find proportional enhancers-gene pairs per library
        lib_df = lib_df.loc[np.sign(lib_df[f"{lib_name}_log2FoldChange_act"])==np.sign(lib_df[f"{lib_name}_log2FoldChange_exp"])]
        lib_save_file = os.path.join(lib_links_save_dir, lib_name, "nearest_sda_sde_table_peaks_proportional.csv")
        os.makedirs(os.path.dirname(lib_save_file), exist_ok=True)
        lib_df.to_csv(lib_save_file, index=False)

