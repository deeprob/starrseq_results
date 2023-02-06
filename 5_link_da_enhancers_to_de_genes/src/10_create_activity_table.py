import argparse
import os
import utils as ut


def create_activity_table(
    depth_dir,
    bam_dir, 
    lib_short,
    enhancerlist,
    store_dir
    ):
    depth_file = os.path.join(depth_dir, lib_short, f"{lib_short}.bed")
    bam_file = os.path.join(bam_dir, lib_short, f"{lib_short}.bam")
    store_file = os.path.join(store_dir, lib_short, "EnhancerActivity.tsv")
    os.makedirs(os.path.dirname(store_file), exist_ok=True)
    # normalize raw counts by RPM and store the annotations learnt by ABC model
    ut.create_activity_table_helper(depth_file, bam_file, enhancerlist, store_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ABC model input preparation")
    parser.add_argument("depth_dir", type=str, help="Dir where the fragment coverage files are stored")
    parser.add_argument("bam_dir", type=str, help="Dir where the filtered bam files are stored")
    parser.add_argument("lib_short", type=str, help="Lib short form")
    parser.add_argument("enhancerlist", type=str, help="Enhancer list file created by ABC model")
    parser.add_argument("store_dir", type=str, help="Output dir where parsed file will be stored")

    cli_args = parser.parse_args()

    create_activity_table(cli_args.depth_dir, cli_args.bam_dir, cli_args.lib_short, cli_args.enhancerlist, cli_args.store_dir)
