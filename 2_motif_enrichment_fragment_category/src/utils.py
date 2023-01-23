import os
import json
from argparse import Namespace
import shutil
import pandas as pd
import pybedtools
import subprocess
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))


###############################
# read meta file; create args #
###############################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)
        
    args = Namespace(
        # from metadata file
        library_prefix = meta_dict[lib_name]["prefix"],
        library_reps = meta_dict[lib_name]["replicates"],
        library_pair= meta_dict[lib_name]["read_pairs"],
        library_umi = meta_dict[lib_name]["umi"],
        library_suffix = meta_dict[lib_name]["suffix"],
        library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["filtered"]
    )

    return args


###################
# filepath parser #
###################

def get_lib_peak_filepath(peak_dir, lib_short, peak_caller="starrpeaker"):
    peak_filename_dict = {"starrpeaker": "peaks.peak.final.bed"}
    peak_filepath = os.path.join(
        peak_dir, lib_short, peak_caller, peak_filename_dict[peak_caller]
        )
    return peak_filepath

def get_lib_dapeak_filepath(peak_dir, lib_short, da_type):
    peak_filepath = os.path.join(
        peak_dir, lib_short, f"{da_type}.bed"
        )
    return peak_filepath

def get_mea_outdir_path(store_dir, analysis, lib_short, da_type, method="homer", preprocess=False):
    suff = "in" if preprocess else "out"
    return os.path.join(store_dir, analysis, lib_short, da_type, f"{method}_{suff}")

def get_background_filepath(preprocess_outdir):
    return os.path.join(preprocess_outdir, f"background.bed")

def get_filebasename(filename):
    return os.path.splitext(os.path.basename(filename))[0]


####################
# homer preprocess #
####################

def create_homer_compatible_bedfile(bed_filepath, processed_filepath):
    df = pd.read_csv(bed_filepath, sep="\t", usecols=[0,1,2], header=None).drop_duplicates(keep="first")
    df[3] = df[0] + ":" + df[1].astype(str) + "-" + df[2].astype(str)
    df[4] = 0
    df[5] = "."
    os.makedirs(os.path.dirname(processed_filepath), exist_ok=True)
    df.to_csv(processed_filepath, sep="\t", header=None, index=None)
    return

def make_windows(in_bed, out_bed, window_size=500, window_stride=50):
    """
    Break the ROIs into fragments of an user defined window size and stride
    """
    window = pybedtools.BedTool().window_maker(b=in_bed ,w=window_size, s=window_stride)
    window_df = window.to_dataframe()
    # get rid of windows which have the same end point
    last_end = None
    rows_to_omit = []
    for i, row in enumerate(window_df.itertuples()):
        if row.end == last_end:
            rows_to_omit.append(i)
        last_end = row.end
    window_df = window_df.loc[~window_df.index.isin(rows_to_omit)]
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    window_df.to_csv(out_bed, sep="\t", header=None, index=None)
    return


def create_background_helper(peak_bed, master_bed, background_path):
    # get background regions 
    background_bed = master_bed - peak_bed
    # save background_windows
    make_windows(background_bed, background_path)
    create_homer_compatible_bedfile(background_path, background_path)
    return


def create_background(master_file, peak_file, background_file):
    master_bed = pybedtools.BedTool(master_file)
    peak_bed = pybedtools.BedTool(peak_file)
    create_background_helper(peak_bed, master_bed, background_file)
    return


#############
# homer mea #
#############

def run_mea_homer(peak_file, reference_genome, background_region, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    logfile = os.path.join(output_dir, "homer.log")

    lf = open(logfile, "w")
    subprocess.run([
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/0_run_mea_homer.sh", 
        peak_file, reference_genome, background_region, output_dir, f"{threads}"
        ], stdout=lf, stderr=lf, check=True)
    lf.close()
    return logfile

###################
# meme preprocess #
###################

def fastafrombed(bed_filepath, genome_filepath, output_filepath):
    ip_file = pybedtools.BedTool(bed_filepath)
    ip_fasta = ip_file.sequence(fi=genome_filepath)
    shutil.move(ip_fasta.seqfn, output_filepath)
    return

def meme_preprocess(preprocess_outdir, peak_file, background_file, reference_genome):
    peak_fasta_file, bg_fasta_file = os.path.join(preprocess_outdir, f"{get_filebasename(peak_file)}.fa"), os.path.join(preprocess_outdir, f"{get_filebasename(background_file)}.fa")
    fastafrombed(peak_file, reference_genome, peak_fasta_file)
    fastafrombed(background_file, reference_genome, bg_fasta_file)
    return peak_fasta_file, bg_fasta_file

############
# meme mea #
############

def run_mea_meme(peak_file, reference_genome, background_region, motif_file, preprocess_outdir, results_outdir):
    peak_fasta_file, bg_fasta_file = meme_preprocess(preprocess_outdir, peak_file, background_region, reference_genome)
    
    os.makedirs(results_outdir, exist_ok=True)
    logfile = os.path.join(results_outdir, "meme.log")

    lf = open(logfile, "w")
    subprocess.run([
        "bash", 
        f"{CURRENT_DIR_PATH}/shell_scripts/0_run_mea_meme.sh", 
        peak_fasta_file, bg_fasta_file, motif_file, results_outdir
        ], stdout=lf, stderr=lf, check=True)
    lf.close()
    return logfile


###########
# mea run #
###########

def run_mea(
    method, 
    lib_peak_file, 
    genome, 
    background_region_filepath,
    motif_file,
    preprocess_outdir, 
    results_outdir, 
    threads
    ):

    if method == "homer":
        run_mea_homer(lib_peak_file, genome, background_region_filepath, results_outdir, threads)
    else:
        # assuming alternate method is meme
        run_mea_meme(lib_peak_file, genome, background_region_filepath, motif_file, preprocess_outdir, results_outdir)
    return


####################
# homer motif scan #
####################

def pwm_scan_homer(homer_roi, homer_motifs, genome, homer_outdir, threads):
    logfile = os.path.join(homer_outdir, "motif_scan.log")
    outfile = os.path.join(homer_outdir, "motif_scan.tsv")
    cmd = [
        "findMotifsGenome.pl", homer_roi, genome, homer_outdir, 
        "-find", homer_motifs, "-p", str(threads)
        ]
    with open(logfile, "w") as lf:
        with open(outfile, "w") as of:
            results = subprocess.run(cmd, stdout=of, stderr=lf)
    return results
