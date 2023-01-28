import argparse
import os
import utils as ut


def create_expression_table(
    root_dir, 
    lib_short, 
    counts_filebase,
    parsed_gtf,
    store_dir
    ):
    counts_file = os.path.join(root_dir, lib_short, counts_filebase)
    store_file = os.path.join(store_dir, lib_short, "expression.tsv")
    os.makedirs(os.path.dirname(store_file), exist_ok=True)
    # normalize raw counts by TPM normalization
    ut.create_expression_table_helper(counts_file, parsed_gtf, store_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ABC model input preparation")
    parser.add_argument("counts_dir", type=str, help="Dir where the rnaseq counts files are stored")
    parser.add_argument("lib_short", type=str, help="Lib short form")
    parser.add_argument("gtf_file", type=str, help="parsed GTF filepath")
    parser.add_argument("store_dir", type=str, help="Output dir where parsed file will be stored")
    parser.add_argument("--counts_filebase", default="counts.tsv", type=str, help="Counts file basename")

    cli_args = parser.parse_args()

    create_expression_table(cli_args.counts_dir, cli_args.lib_short, cli_args.counts_filebase, cli_args.gtf_file, cli_args.store_dir)
