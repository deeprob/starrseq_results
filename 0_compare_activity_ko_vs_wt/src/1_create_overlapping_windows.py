import os
import pandas as pd
import argparse
import utils as ut


def main(
    master_list_dir,
    store_dir,
    lib_short,
    faster
    ):
    # get the master file filtered with at least 50 reads per replicate
    master_file = os.path.join(master_list_dir, lib_short, "master_filtered.bed")
    # window out file
    store_file = os.path.join(store_dir, "window", lib_short, "master_filtered_windows.bed")
    os.makedirs(os.path.dirname(store_file), exist_ok=True)
    # create windows
    if faster:
        ut.make_windows_faster(master_file, store_file)
    else:
        ut.make_windows(master_file, store_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib_name", type=str, help="The library name as given in the meta file")
    parser.add_argument("master_file_dir", type=str, help="Dir where the master list is stored")
    parser.add_argument("store_dir", type=str, help="Save file dir")
    parser.add_argument("-f", "--faster", action="store_true", help="Use if the roi file is big")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib_name)

    main(
        cli_args.master_file_dir,
        cli_args.store_dir,
        lib_args.library_short,
        cli_args.faster        
    )