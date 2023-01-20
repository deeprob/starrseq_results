import os
import argparse
import utils as ut


def main(
    filtered_bam_dir,
    store_dir,
    window_file,
    lib_prefix,
    lib_replicates,
    lib_short,
    ):
    # get the bam file
    replicate_bam_files = [os.path.join(filtered_bam_dir, lib_short, f"{lib_prefix}_{lib_rep}.bam") for lib_rep in lib_replicates.split()]
    merged_bam_file = os.path.join(filtered_bam_dir, lib_short, f"{lib_prefix}.bam")
    # depth out file
    replicate_depth_files = [os.path.join(store_dir, "window_depth", lib_short, f"{lib_prefix}_{lib_rep}.bed") for lib_rep in lib_replicates.split()]
    merged_depth_file = os.path.join(store_dir, "window_depth", lib_short, f"{lib_prefix}.bed")
    os.makedirs(os.path.dirname(merged_depth_file), exist_ok=True)
    # create windows
    pool_iter = [(bf, window_file, sf) for bf,sf in zip(replicate_bam_files+[merged_bam_file], replicate_depth_files + [merged_depth_file])]
    ut.run_multiargs_pool_job(ut.get_roi_depth, pool_iter)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib_name", type=str, help="The library name as given in the meta file")
    parser.add_argument("filtered_bam_dir", type=str, help="Dir where the filtered bam files are stored")
    parser.add_argument("store_dir", type=str, help="Save files dir")
    parser.add_argument("roi_window_file", type=str, help="Path to the roi window file whose library coverage will be estimated")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib_name)

    main(
        cli_args.filtered_bam_dir,
        cli_args.store_dir,
        cli_args.roi_window_file,
        lib_args.library_prefix, 
        lib_args.library_reps,
        lib_args.library_short        
    )
