import argparse
import os
import utils as ut


def parse_gtf(gtf_file, save_dir):
    parsed_gtf_file = os.path.join(save_dir, "genome_annot", "parsed_gtf.tsv")
    os.makedirs(os.path.dirname(parsed_gtf_file), exist_ok=True)
    ut.parse_gtf_file(gtf_file, parsed_gtf_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ABC model input preparation")
    parser.add_argument("gtf_file", type=str, help="GTF filepath to parse")
    parser.add_argument("store_dir", type=str, help="Output dir where parsed file will be stored")

    cli_args = parser.parse_args()

    parse_gtf(cli_args.gtf_file, cli_args.store_dir)
