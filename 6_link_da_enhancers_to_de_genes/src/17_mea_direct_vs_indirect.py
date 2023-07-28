import os
import numpy as np
import pandas as pd


def create_homer_bed(df, libn, direct_save_file, indirect_save_file, direct_loss_save_file, indirect_gained_save_file):
    if len(df)==0:
        return
    # drop duplicates
    df = df.fillna("NA").drop_duplicates(keep="first")
    # split the merged coordinates to bed format
    df = pd.concat((df, df.chrom_coord.str.split("_", expand=True).rename(columns={0: "chrom", 1: "start", 2: "end"})), axis=1)
    # add strand info
    df["strand"] = "."
    df = df.loc[:, ["chrom", "start", "end", "chrom_coord", f"{libn}_act", "strand", f"{libn}_log2FoldChange_act", f"gene_name", f"{libn}_log2FoldChange_exp", f"{libn}_padj_act", f"{libn}_padj_exp", "CC_peak", f"{libn}_peak"]]
    direct_df = df.loc[(df[f"{libn}_log2FoldChange_act"]<0)&(df[f"{libn}_log2FoldChange_exp"]<0)]
    direct_loss_df = df.loc[(df[f"CC_peak"]==1)&(df[f"{libn}_peak"]==0)&(df[f"{libn}_log2FoldChange_act"]<0)&(df[f"{libn}_log2FoldChange_exp"]<0)]
    indirect_df = df.loc[(df[f"{libn}_log2FoldChange_act"]>0)&(df[f"{libn}_log2FoldChange_exp"]>0)]
    indirect_gained_df = df.loc[(df[f"CC_peak"]==0)&(df[f"{libn}_peak"]==1)&(df[f"{libn}_log2FoldChange_act"]>0)&(df[f"{libn}_log2FoldChange_exp"]>0)]
    # save to bed format
    direct_df.to_csv(direct_save_file, sep="\t", header=False, index=False)
    direct_loss_df.to_csv(direct_loss_save_file, sep="\t", header=False, index=False)
    indirect_df.to_csv(indirect_save_file, sep="\t", header=False, index=False)
    indirect_gained_df.to_csv(indirect_gained_save_file, sep="\t", header=False, index=False)
    return


if __name__=="__main__":
    sda_sde_table_dir = "/data5/deepro/starrseq/papers/results/6_link_da_enhancers_to_de_genes/data/da_de_peaks"
    libraries = ["ATF2", "CTCF", "FOXA1", "LEF1", "SCRT1", "TCF7L2", "16P12_1"]
    for lib_name in libraries:
        abc_lib_file = os.path.join(sda_sde_table_dir, lib_name, "abc_sda_sde_table_peaks_fragments.csv")
        abc_lib_df = pd.read_csv(abc_lib_file)
        nearest_lib_file = os.path.join(sda_sde_table_dir, lib_name, "nearest_sda_sde_table_peaks_fragments.csv")
        nearest_lib_df = pd.read_csv(nearest_lib_file)
        lib_df = pd.concat((abc_lib_df, nearest_lib_df)).drop_duplicates(subset=["chrom_coord", "gene_name"])
        dsave_file = os.path.join(sda_sde_table_dir, lib_name, "direct.bed")
        idsave_file = os.path.join(sda_sde_table_dir, lib_name, "indirect.bed")
        dlsave_file = os.path.join(sda_sde_table_dir, lib_name, "direct_loss.bed")
        idgsave_file = os.path.join(sda_sde_table_dir, lib_name, "indirect_gained.bed")
        create_homer_bed(lib_df, lib_name, dsave_file, idsave_file, dlsave_file, idgsave_file)
