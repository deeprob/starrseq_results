import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns


def create_activity_expression_df(meta_activity_df, meta_expression_df, target_gene_col, libraries, act_column_suff, exp_column_suff):
    activity_df = meta_activity_df.loc[~meta_activity_df[target_gene_col].isna()]
    activity_df[target_gene_col] = activity_df[target_gene_col].str.split("|")
    activity_df = activity_df.explode(target_gene_col)
    activity_df = activity_df.merge(meta_expression_df, left_on=target_gene_col, right_on="gene_id", suffixes=('_act', '_exp'))
    return activity_df

def get_per_gene_corr(ser, ko_lines):
    act_val = ser.loc[[f"{ko}_act" for ko in ko_lines]].values
    exp_val = ser.loc[[f"{ko}_exp" for ko in ko_lines]].values
    return pearsonr(act_val, exp_val)[0]

def check_sde(ser, meta_exp, libs):
    coord = ser.chrom_coord
    genes = np.array(ser.nearest_gene.split("|"))
    sdes = meta_exp.loc[meta_exp.gene_id.isin(genes), [f"{lib}_padj" for lib in libs]].values
    genes = genes.reshape(len(genes),1).repeat(len(libs), 1)
    sde_genes = np.where(sdes<0.01, genes, "")
    data_dict = dict()
    for i,row in enumerate(sde_genes.T):
        data_dict[f"{libs[i]}_sde"] = "|".join(row).strip("|")
    return pd.Series(data_dict)


if __name__ == "__main__":
    meta_activity_df = pd.read_csv("/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/meta_enhancer_gene.csv")
    meta_expression_df = pd.read_csv("/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/meta_exp.csv")
    libraries = ["CC", "ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    abc_save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/abc_da_de_table.csv"
    nearest_save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/nearest_da_de_table.csv"
    sde_save_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/activity_vs_expression_corr/nearest_da_sde_table.csv"

    abc_activity_df = create_activity_expression_df(meta_activity_df, meta_expression_df, "abc_gene", libraries[1:], "log2FoldChange", "log2FoldChange")
    nearest_activity_df = create_activity_expression_df(meta_activity_df, meta_expression_df, "nearest_gene", libraries[1:], "log2FoldChange", "log2FoldChange")
    abc_activity_df["per_gene_corr"] = abc_activity_df.apply(get_per_gene_corr, args=(libraries[1:],), axis=1)
    nearest_activity_df["per_gene_corr"] = nearest_activity_df.apply(get_per_gene_corr, args=(libraries[1:],), axis=1)
    abc_activity_df.sort_values("per_gene_corr", ascending=False).to_csv(abc_save_file)
    nearest_activity_df.sort_values("per_gene_corr", ascending=False).to_csv(nearest_save_file)

    
    meta_activity_df = pd.concat([meta_activity_df, meta_activity_df.apply(check_sde,  args=(meta_expression_df, libraries[1:]), axis=1)], axis=1)
    meta_activity_df.to_csv(sde_save_file)
