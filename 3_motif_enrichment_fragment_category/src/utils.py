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
        roi_file = meta_dict["roi"]["sorted"]
    )

    return args


#############
# homer mea #
#############

def run_mea_homer(peak_file, reference_genome, background_region, motif_file, output_dir, threads):
    os.makedirs(output_dir, exist_ok=True)
    logfile = os.path.join(output_dir, "homer.log")

    lf = open(logfile, "w")
    subprocess.run([
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/0_run_mea_homer.sh", 
        peak_file, reference_genome, background_region, output_dir, f"{threads}", motif_file
        ], stdout=lf, stderr=lf, check=True)
    lf.close()
    return logfile

###################
# meme preprocess #
###################

def get_filebasename(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def fastafrombed(bed_filepath, genome_filepath, output_filepath):
    ip_file = pybedtools.BedTool(bed_filepath)
    ip_fasta = ip_file.sequence(fi=genome_filepath)
    shutil.move(ip_fasta.seqfn, output_filepath)
    return

def meme_preprocess(outdir, peak_file, background_file, reference_genome):
    peak_fasta_file, bg_fasta_file = os.path.join(outdir, f"{get_filebasename(peak_file)}.fa"), os.path.join(outdir, f"{get_filebasename(background_file)}.fa")
    fastafrombed(peak_file, reference_genome, peak_fasta_file)
    fastafrombed(background_file, reference_genome, bg_fasta_file)
    return peak_fasta_file, bg_fasta_file

############
# meme mea #
############

def run_mea_meme(peak_file, reference_genome, background_region, motif_file, outdir):
    os.makedirs(outdir, exist_ok=True)
    peak_fasta_file, bg_fasta_file = meme_preprocess(outdir, peak_file, background_region, reference_genome)
    
    logfile = os.path.join(outdir, "meme.log")

    lf = open(logfile, "w")
    subprocess.run([
        "bash", 
        f"{CURRENT_DIR_PATH}/shell_scripts/0_run_mea_meme.sh", 
        peak_fasta_file, bg_fasta_file, motif_file, outdir
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
    outdir,
    threads
    ):

    if method == "homer":
        run_mea_homer(lib_peak_file, genome, background_region_filepath, motif_file, outdir, threads)
    else:
        # assuming alternate method is meme
        run_mea_meme(lib_peak_file, genome, background_region_filepath, motif_file, outdir)
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
