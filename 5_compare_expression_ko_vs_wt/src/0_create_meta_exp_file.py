import os
import pandas as pd


def get_tpm_normalized_gene_expression(filename, lib):
    df = pd.read_csv(filename, sep="\t", index_col=[0, 1], skipfooter=5, engine="python")
    df = (df*1e6/df.sum()).mean(axis=1)
    df.name = lib
    return df.reset_index()

def get_diff_exp(filename, lib):
    df = pd.read_csv(filename, usecols=["log2FoldChange", "pvalue", "padj"])
    df.columns = [f"{lib}_{c}" for c in df.columns]
    return df.reset_index().rename(columns={"index":"gene_id"})


def get_meta_gene_exp_df(exp_dir, de_dir, libraries, save_file):
    df_list = []
    for ko in libraries:
        exp_file = os.path.join(exp_dir, ko, "counts.tsv")
        exp_df = get_tpm_normalized_gene_expression(exp_file, ko)
        diff_exp_file = os.path.join(de_dir, f"{ko}vsCC", "de_results.csv")
        if ko != "CC":
            diff_exp_df = get_diff_exp(diff_exp_file, ko)
            exp_df = exp_df.merge(diff_exp_df, left_on="gene_id", right_on="gene_id")
        exp_df = exp_df.set_index(["gene_id", "gene_name"])
        df_list.append(exp_df)
    exp_df = pd.concat(df_list, axis=1)
    exp_df.to_csv(save_file, index=True)
    return 


if __name__ == "__main__":
    exp_dir = "/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/results/count"
    de_dir = "/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/results/de"
    save_file = "/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/meta_exp.csv"
    libraries = libraries = ["CC", "ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    get_meta_gene_exp_df(exp_dir, de_dir, libraries, save_file)

