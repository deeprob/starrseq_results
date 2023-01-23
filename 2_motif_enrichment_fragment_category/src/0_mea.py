import argparse
import utils as ut


def main(
    master_file,
    lib_short,
    genome,
    peak_store_dir,
    mea_out_dir,
    method,
    motif_file,
    diff_activity_type,
    diff_activity_analysis_type,
    threads,
    ):

    # determine analysis type
    analysis = "lib_peaks"
    if diff_activity_type:
        analysis = diff_activity_analysis_type

    # get peak filepath
    lib_peak_file = ut.get_lib_peak_filepath(peak_store_dir, lib_short) if not diff_activity_type else ut.get_lib_dapeak_filepath(peak_store_dir, lib_short, diff_activity_type)
    
    # get outdir
    preprocess_outdir = ut.get_mea_outdir_path(mea_out_dir, analysis, lib_short, diff_activity_type, method=method, preprocess=True)
    results_outdir = ut.get_mea_outdir_path(mea_out_dir, analysis, lib_short, diff_activity_type, method=method, preprocess=False)

    # create background regions
    background_region_filepath = ut.get_background_filepath(preprocess_outdir)
    ut.create_background(master_file, lib_peak_file, background_region_filepath)

    # conduct mea analysis with homer or meme

    ut.run_mea(
        method, 
        lib_peak_file, 
        genome, 
        background_region_filepath,
        motif_file,
        preprocess_outdir, 
        results_outdir, 
        threads
        )

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq MEA analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib", type=str, help="library name as given in the meta file")
    parser.add_argument("peak_dir", type=str, help="Dir where peak files are stored")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")
    parser.add_argument("--method", type=str, default="homer", help="MEA method - homer/meme")
    parser.add_argument("--motif_file", type=str, default="", help="known motifs to be used with MEA method")
    parser.add_argument("--diff_activity_type", type=str, default="", help="type of differential enhancer activity, use this argument to compare induced,repressed and constitutive peaks between lib1 and control")
    parser.add_argument("--diff_activity_analysis_type", type=str, default="diff_peaks", help="type of method used to measure differential enhancer activity, use this argument to compare differential peaks called using deseq2 or normal bedtools intersection")
    parser.add_argument("--threads", type=int, default=64, help="Number of cores to use")

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)

    main(
        lib_args.roi_file,
        lib_args.library_short,
        lib_args.reference_genome,
        cli_args.peak_dir,
        cli_args.store_dir,
        cli_args.method,
        cli_args.motif_file,
        cli_args.diff_activity_type,
        cli_args.diff_activity_analysis_type,
        cli_args.threads,
    )
