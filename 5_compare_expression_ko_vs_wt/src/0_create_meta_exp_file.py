import os
import pandas as pd


def get_tpm_normalized_gene_expression(counts_file, gtf_file, lib):
    expr_df = pd.read_csv(counts_file, sep="\t", engine="python", skipfooter=5).set_index(["gene_id", "gene_name"])
    expr_df["mean_raw_counts"] = expr_df.mean(axis=1)
    gtf_df = pd.read_csv(gtf_file, sep="\t", header=None, names=["chrm", "start", "end", "gene_id", "gene_name", "strand"]).set_index(["gene_id", "gene_name"])
    expr_df = expr_df.merge(gtf_df, left_index=True, right_index=True)
    assert len(expr_df) == len(gtf_df)
    expr_df["length"] = expr_df.end - expr_df.start
    expr_df["expression"] = (expr_df.mean_raw_counts*1e3)/expr_df.length
    expr_df["expression"] = (expr_df["expression"]*1e6)/expr_df["expression"].sum()
    expr_df = expr_df.reset_index()
    expr_df = expr_df.rename(columns={"expression": lib})
    return expr_df.loc[:, ["gene_id", "gene_name", lib]]


def get_diff_exp(filename, lib):
    df = pd.read_csv(filename, usecols=["log2FoldChange", "pvalue", "padj"])
    df.columns = [f"{lib}_{c}" for c in df.columns]
    return df.reset_index().rename(columns={"index":"gene_id"})


def get_meta_gene_exp_df(exp_dir, de_dir, libraries, geno_annot_file, save_file):
    df_list = []
    for ko in libraries:
        exp_file = os.path.join(exp_dir, ko, "counts.tsv")
        exp_df = get_tpm_normalized_gene_expression(exp_file, geno_annot_file, ko)
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
    geno_annot_file = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/genome_annot/parsed_gtf.tsv"
    save_file = "/data5/deepro/starrseq/papers/results/5_compare_expression_ko_vs_wt/data/meta_exp.csv"
    libraries = libraries = ["CC", "ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    get_meta_gene_exp_df(exp_dir, de_dir, libraries, geno_annot_file, save_file)

