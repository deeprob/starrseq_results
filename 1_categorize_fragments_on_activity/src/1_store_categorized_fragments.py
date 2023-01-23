import os
import argparse
import utils as ut
import pandas as pd


def save_categorized_fragments(meta_activity_map_file, all_lib_names, store_dir):
    meta_df = pd.read_csv(meta_activity_map_file)
    for libn in all_lib_names:
        if libn == "control":
            ut.save_peaks_notpeaks(meta_df, libn, store_dir, True)        
        else:
            ut.save_peaks_notpeaks(meta_df, libn, store_dir, False)
            ut.save_diff_act(meta_df, libn, store_dir)
    ut.save_always_active_inactive(meta_df, all_lib_names, store_dir)
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
