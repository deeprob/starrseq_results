import os
import argparse
import utils as ut


def main(
    lib_short,
    lib_prefix,
    lib_reps, 
    cc_short,
    cc_prefix,
    cc_reps,
    depth_dir,
    store_dir
    ):
    store_dir = os.path.join(store_dir, "diff_activity")
    lib_depth_files = ut.get_rep_depth_files(depth_dir, lib_short, lib_prefix, lib_reps)
    cc_depth_files = ut.get_rep_depth_files(depth_dir, cc_short, cc_prefix, cc_reps)
    ut.get_diff_active_fragments_deseq(lib_depth_files, cc_depth_files, lib_short, cc_short, store_dir)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq differentially active regions")
    parser.add_argument("meta_file", type=str, help="The file path to the meta file")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("control_lib", type=str, default="", help="library name of the control line as given in the meta file, compare induced,repressed and non-responsive peaks between ko and control")
    parser.add_argument("depth_dir", type=str, help="Dir where windowed depth files are stored")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    cc_args = ut.create_args(cli_args.meta_file, cli_args.control_lib)

    main(
        lib_args.library_short,
        lib_args.library_prefix,
        lib_args.library_reps,
        cc_args.library_short,
        cc_args.library_prefix,
        cc_args.library_reps,
        cli_args.depth_dir,
        cli_args.store_dir
    )
