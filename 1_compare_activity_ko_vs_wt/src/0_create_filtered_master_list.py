import os
import pandas as pd
import argparse
import utils as ut


def read_and_extract_coverage(cov_bed):
    df = pd.read_csv(cov_bed, sep="\t", header=None)
    roi_depth = df.iloc[:, [0,1,2,-4]].set_index([0,1,2])
    return roi_depth

def main(
    cov_dir,
    store_dir,
    lib_prefix,
    lib_replicates,
    lib_short,
    ):
    lib_depth_beds = ut.get_lib_depth_beds_filepaths(cov_dir, lib_short, lib_prefix, lib_replicates)
    df_cov_beds = pd.concat(list(map(read_and_extract_coverage, lib_depth_beds)), axis=1)
    # only keep those regions with greater than 50 reads across all replicates
    df_cov_beds = df_cov_beds.loc[(df_cov_beds>50).all(axis=1)]
    store_file = os.path.join(store_dir, "master", lib_short, "master_filtered.bed")
    os.makedirs(os.path.dirname(store_file), exist_ok=True)
    df_cov_beds.to_csv(store_file, sep="\t", header=None, index=True)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib_name", type=str, help="The library name as given in the meta file")
    parser.add_argument("cov_dir", type=str, help="Coverage files dir")
    parser.add_argument("store_dir", type=str, help="Save files dir")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib_name)

    main(
        cli_args.cov_dir,
        cli_args.store_dir,
        lib_args.library_prefix, 
        lib_args.library_reps,
        lib_args.library_short        
    )