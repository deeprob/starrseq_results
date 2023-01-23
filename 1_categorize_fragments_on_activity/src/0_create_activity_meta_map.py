import os
import argparse
import utils as ut


def save_meta_activity_map(meta_file, all_lib_names, fragment_depth_dir, filtered_bam_dir, peak_dir, diff_act_dir, store_dir):
    meta_df = ut.get_meta_activity_map(meta_file, all_lib_names, fragment_depth_dir, filtered_bam_dir, peak_dir, diff_act_dir)
    save_file = os.path.join(store_dir, "meta_activity_map.csv")
    meta_df.to_csv(save_file, index=False)
    return


if  __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("window_cov_dir", type=str, help="Fragment coverage files dir")
    parser.add_argument("bam_dir", type=str, help="Filtered bam dir for the libraries")
    parser.add_argument("peak_dir", type=str, help="Starrpeaker peak dir for the libraries")
    parser.add_argument("diff_act_dir", type=str, help="Deseq2 diff act dir for the libraries")
    parser.add_argument("store_dir", type=str, help="Save files dir")
    parser.add_argument("-l", "--lib_names", type=str, help="The library names as given in the meta file", nargs="+")

    cli_args = parser.parse_args()
    
    save_meta_activity_map(
        cli_args.meta_file, 
        cli_args.lib_names, 
        cli_args.window_cov_dir, 
        cli_args.bam_dir, 
        cli_args.peak_dir, 
        cli_args.diff_act_dir, 
        cli_args.store_dir
    )
