import os
import argparse
import utils as ut


def main(
    lib_short,
    peak_file,
    background_file,
    genome,
    motif_file,
    store_dir,
    comparison,
    method,
    threads,
    ):

    store_dir = os.path.join(store_dir, lib_short, method, comparison)

    # conduct mea analysis with homer or meme
    ut.run_mea(
        method, 
        peak_file, 
        genome, 
        background_file,
        motif_file, 
        store_dir,
        threads
        )

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq MEA analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("peak_file", type=str, help="MEA peak file path")
    parser.add_argument("control_file", type=str, help="MEA background file path")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")
    parser.add_argument("comparison", type=str, help="Types being compared")
    parser.add_argument("--method", type=str, default="homer", help="MEA method - homer/meme")
    parser.add_argument("--motif_file", type=str, default="", help="known motifs to be used with MEA method")
    parser.add_argument("--threads", type=int, default=64, help="Number of cores to use")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)

    main(
        lib_args.library_short,
        cli_args.peak_file,
        cli_args.control_file,
        lib_args.reference_genome,
        cli_args.motif_file,
        cli_args.store_dir,
        cli_args.comparison,
        cli_args.method,
        cli_args.threads,
    )
