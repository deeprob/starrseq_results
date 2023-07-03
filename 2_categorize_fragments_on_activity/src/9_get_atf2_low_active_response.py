import pandas as pd


def get_atf2_anomaly(meta_activity_file, save_file):
    meta_df = pd.read_csv(meta_activity_file, low_memory=False)
    cc_decile_ranks, decile_labels = pd.qcut(meta_df.CC, q=10, labels=False, retbins=True)
    meta_df["cc_decile"] = cc_decile_ranks
    meta_df[f"ATF2_response"] = meta_df[f"ATF2_padj"]<0.01 
    atf2_df = meta_df.loc[(meta_df.ATF2_response==1)&(meta_df.cc_decile<5)]
    atf2_df = atf2_df.merge(atf2_df.chrom_coord.str.split("_", expand=True), left_index=True, right_index=True).rename(columns={0:"chrm", 1:"start", 2:"end"})
    atf2_df["strand"] = "."
    atf2_df = atf2_df.loc[:, ["chrm", "start", "end", "chrom_coord", f"ATF2", "strand", f"ATF2_log2FoldChange", f"ATF2_baseMean", f"ATF2_lfcSE", f"ATF2_pvalue", f"ATF2_padj"]]
    atf2_df.to_csv(save_file, index=False, sep="\t", header=False)
    return

if __name__ == "__main__":
    meta_activity_file = "/data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/data/meta_activity_map.csv"
    save_file = "/data5/deepro/starrseq/papers/results/2_categorize_fragments_on_activity/data/ATF2/responsive_low_activity.bed"
    get_atf2_anomaly(meta_activity_file, save_file)
