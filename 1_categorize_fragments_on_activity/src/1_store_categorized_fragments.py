import os
import argparse
import utils as ut
import pandas as pd


def save_categorized_fragments(meta_activity_map_file, all_lib_names, store_dir):
    meta_df = pd.read_csv(meta_activity_map_file)
    for libn in all_lib_names:
        
        meta_df[f"{libn}_padj"] = meta_df[f"{libn}_padj"].fillna(1.)
        # induced fragments :: ko rpp activity > 0; log2FoldChange>0; padj<0.01 
        induced = meta_df.loc[(meta_df[libn]>0)&(meta_df[f"{libn}_log2FoldChange"]>0)&(meta_df[f"{libn}_padj"]<0.01)]
        ut.save_diff_act_file(induced, store_dir, libn, "induced")
        # repressed fragments :: control rpp activity > 0; log2FoldChange<0; padj<0.01 
        repressed = meta_df.loc[(meta_df["control"]>0)&(meta_df[f"{libn}_log2FoldChange"]<0)&(meta_df[f"{libn}_padj"]<0.01)]
        ut.save_diff_act_file(repressed, store_dir, libn, "repressed")
        # non-responsive fragments :: control or ko rpp activity > 0; padj>0.01 
        nonresponsive = meta_df.loc[((meta_df["control"]>0)|(meta_df[libn]>0))&(meta_df[f"{libn}_padj"]>0.01)]
        ut.save_diff_act_file(nonresponsive, store_dir, libn, "nonresponsive")
    return


if  __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_act_file", type=str, help="Filepath to the meta activity map")
    parser.add_argument("store_dir", type=str, help="Save files dir")
    parser.add_argument("-l", "--lib_names", type=str, help="The library names as given in the meta file", nargs="+")

    cli_args = parser.parse_args()
    
    save_categorized_fragments(
        cli_args.meta_act_file,
        cli_args.lib_names, 
        cli_args.store_dir
    )
