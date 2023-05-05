import os
import numpy as np
import pandas as pd


def get_fragment_categories(df, lib):
    # filter to keep only CC and lib info
    relevant_columns = ["chrom_coord", "CC", "CC_peak"] + [col for col in df.columns if col.startswith(lib)]
    df = df.loc[:, relevant_columns]

    unresponsive_cond = f"(`{lib}_padj` > 0.01)"
    induced_cond = f"(`{lib}_padj` < 0.01) & (`{lib}_log2FoldChange` > 0)"
    repressed_cond = f"(`{lib}_padj` < 0.01) & (`{lib}_log2FoldChange` < 0)"
    active_induced_cond = "(`CC_peak` == 1) & " f"(`{lib}_padj` < 0.01) & (`{lib}_peak` == 1) & (`{lib}_log2FoldChange` > 0)"
    gained_cond = "(`CC_peak` == 0) & " + f"(`{lib}_padj` < 0.01) & (`{lib}_peak` == 1) & (`{lib}_log2FoldChange` > 0)" 
    lost_cond = "(`CC_peak` == 1) & " + f"(`{lib}_padj` < 0.01) & (`{lib}_peak` == 0) & (`{lib}_log2FoldChange` < 0)"
    always_active_cond =  "(`CC_peak` == 1) & " +  f"(`{lib}_peak` == 1)"
    always_inactive_cond =  "(`CC` < -1) & " +  f"(`{lib}` < -1)" 

    unresponsive_fragments = df.query(unresponsive_cond)
    induced_fragments = df.query(induced_cond)
    repressed_fragments = df.query(repressed_cond)
    active_induced_fragments = df.query(active_induced_cond)
    always_active_fragments = df.query(always_active_cond)
    always_inactive_fragments = df.query(always_inactive_cond)
    gained_fragments = df.query(gained_cond)
    lost_fragments = df.query(lost_cond)

    return unresponsive_fragments, induced_fragments, repressed_fragments, always_active_fragments, always_inactive_fragments, gained_fragments, lost_fragments

def parse_deseqres_for_volcano_plot(df, lib):
    df = df.loc[:, ["gene_id", "gene_name"] + [c for c in df.columns if c.startswith(lib)]]
    df = df.rename(columns={f"{lib}_log2FoldChange": "log2FoldChange", f"{lib}_padj": "padj"})
    # drop rows with na values
    df = df.dropna()
    # convert all 0 padj values to a tenth of the min padj value thats greater than 0
    df.loc[df.padj==0, "padj"] = min(df.loc[df.padj>0].padj)/10
    # create neglog10 val
    df["neglog10padj"] = -np.log10(df.padj)
    # create hue columns
    lfc_thresh = 0.5
    pv_thresh = 0.01
    df["hue"] = "Not Significant"
    df.loc[(df.log2FoldChange>lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Up"
    df.loc[(df.log2FoldChange<-lfc_thresh)&(df.padj<pv_thresh), "hue"] = "Significant Down"
    return df

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
        "CC": max, # max CC activity
        lib: max, # max activity value for each fragment
        f"{lib}_pvalue": min, # min significance value for each fragment
        f"{lib}_padj": min, # min significance value for each fragment
        f"{lib}_sde": lambda x: list(x)[0],
        "CC_peak": lambda x: list(x)[0], 
        f"{lib}_peak": lambda x: list(x)[0], 
        f"{lib}_log2FoldChange": lambda x: max(x) if all(i>0 for i in x) else min(x),
        }
        
    grouped_dfs = []
    for group, df in merge_df.groupby("chr"): # https://stackoverflow.com/questions/46732760/merge-rows-pandas-dataframe-based-on-condition
        grouped_df = df.groupby(((df.start  - df.end.shift(1)) > 0).cumsum()).agg(dict_function)
        grouped_df["chr"] = group
        grouped_dfs.append(grouped_df)

    merged_df = pd.concat(grouped_dfs)
    return merged_df.drop(columns="chr")

def add_expression_info(fragment_df, volcano_df, lib):
    sde_col = f"{lib}_sde"
    fragment_df = fragment_df.loc[~fragment_df[sde_col].isna()]
    fragment_df = merge_multi_fragment_into_one(fragment_df, lib)
    fragment_df[sde_col] = fragment_df[sde_col].str.split("|")
    fragment_df = fragment_df.explode(sde_col)
    fragment_df = fragment_df.merge(volcano_df.loc[:, ["gene_id", "gene_name", "log2FoldChange", "neglog10padj"]], left_on=sde_col, right_on="gene_id")
    return fragment_df.sort_values(["neglog10padj", f"{lib}_padj"], ascending=False)


def add_library_info(meta_activity_df, meta_expression_df, lib, save_dir):
    unresponsive_fragments, induced_fragments, repressed_fragments, always_active_fragments, always_inactive_fragments, gained_fragments, lost_fragments = get_fragment_categories(meta_activity_df, lib)
    volcano_df = parse_deseqres_for_volcano_plot(meta_expression_df, lib)
    save_dir = os.path.join(save_dir, lib)
    os.makedirs(save_dir, exist_ok=True)
    volcano_file = os.path.join(save_dir, "volcano.csv")
    volcano_df.to_csv(volcano_file)

    # add Enhancer to Gene targets
    for fragment_df, fragment_name in zip([unresponsive_fragments, induced_fragments, repressed_fragments, gained_fragments, lost_fragments], ["unresponsive", "induced", "repressed", "gained", "lost"]):
        if len(fragment_df)>0:
            save_file = os.path.join(save_dir, f"{fragment_name}.csv")
            fragment_df.to_csv(save_file, index=False)
            tg_df = add_expression_info(fragment_df, volcano_df, lib)
            save_file = os.path.join(save_dir, f"{fragment_name}_merged.csv")
            tg_df.to_csv(save_file, index=False)
    return


if __name__ == "__main__":
    meta_activity_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/nearest_da_sde_table.csv"
    meta_expression_file = "/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/meta_exp.csv"
    save_dir = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_enhancers_to_de_genes_links"
    libs = ["ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]

    meta_activity_df = pd.read_csv(meta_activity_file, low_memory=False)
    meta_expression_df = pd.read_csv(meta_expression_file)

    for lib in libs:
        add_library_info(meta_activity_df, meta_expression_df, lib, save_dir)
