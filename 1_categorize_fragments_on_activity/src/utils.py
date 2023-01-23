import os
import numpy as np
import pandas as pd
import pybedtools
import pysam
import json
from argparse import Namespace
import itertools
import functools
import multiprocessing as mp

# TODO: take tmp dir from arguments
pybedtools.helpers.set_tempdir("/data5/deepro/tmp/")


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


##############
# fileparser #
##############

def get_diff_act_store_file(store_dir, libname, activity_type, fileext):
    return os.path.join(store_dir, libname, f"{activity_type}.{fileext}")


#####################
# meta activity map #
#####################

def read_and_extract_coverage(cov_bed):
    df = pd.read_csv(cov_bed, sep="\t", header=None, index_col=[0,1,2])
    roi_depth = df.iloc[:, -4]
    return roi_depth

def get_bam_file_reads(bam_file):
    reads = pysam.view("-c", bam_file)
    return int(reads.strip())

def get_bed_annotations_helper(read_df, peak_file):
    # read the dataframe as a bed file
    roi_bed = pybedtools.BedTool.from_dataframe(read_df.chrom_coord.str.split("_", expand=True))
    # read the peak bed file
    peak_bed = pybedtools.BedTool(peak_file)
    # find intersections between roi and peak bed, intersection must be at least 95%
    intersect_df = roi_bed.intersect(peak_bed, c=True, f=0.95).to_dataframe(disable_auto_names=True, header=None)
    return (intersect_df.iloc[:, -1]>0).astype(int)

def get_bed_annotations(read_df, bed_files, libraries, bed_type):
    bed_iter = [(read_df, pf) for pf in bed_files]
    intersect_dfs = list(itertools.starmap(get_bed_annotations_helper, bed_iter))
    intersect_df_meta = pd.concat(intersect_dfs, axis=1)
    intersect_df_meta.columns = [f"{l}_{bed_type}" for l in libraries]
    read_df_with_annots = pd.concat((read_df, intersect_df_meta), axis=1)
    return read_df_with_annots

def read_diff_output(diff_act_filename, libname):
    df = pd.read_csv(diff_act_filename)
    df.columns = [f"{libname}_{c}" for c in df.columns]
    return df

def get_meta_activity_map(meta_file, all_lib_names, fragment_depth_dir, filtered_bam_dir, peak_dir, diff_act_dir):
    # get required file paths
    lib_depth_files = []
    bam_files = []
    peak_files = []
    peak_file_libname = []
    diff_act_files = []
    diff_act_libname = []
    for lib_name in all_lib_names:
        lib_args = create_args(meta_file, lib_name)
        lib_depth_file = os.path.join(fragment_depth_dir, lib_args.library_short, f"{lib_args.library_prefix}.bed")
        lib_depth_files.append(lib_depth_file)
        bam_file = os.path.join(filtered_bam_dir, lib_args.library_short,  f"{lib_args.library_prefix}.bam")
        bam_files.append(bam_file)
        peak_file = os.path.join(peak_dir, lib_args.library_short,  "starrpeaker", "peaks.peak.bed")
        diff_act_file = os.path.join(diff_act_dir, lib_args.library_short, "deseq_out.csv")
        if lib_name == "input":
            pass
        elif lib_name == "control":
            peak_files.append(peak_file)
            peak_file_libname.append(lib_name)
        else:
            peak_files.append(peak_file)
            peak_file_libname.append(lib_name)
            diff_act_files.append(diff_act_file)
            diff_act_libname.append(lib_name)
    # create the rpm normalized rpp dataframe
    read_df = pd.concat(list(map(read_and_extract_coverage, lib_depth_files)), axis=1)
    read_df.index = list(map(lambda x: "_".join(list(map(str, x))) , read_df.index))
    read_df.columns = all_lib_names
    read_df = read_df.loc[read_df.input>150]
    bam_file_reads = run_singleargs_pool_job(get_bam_file_reads, bam_files) # list(map(get_bam_file_reads, bam_files))
    read_rpm_df = read_df.divide(bam_file_reads, axis=1)
    read_rpm_df = read_rpm_df * 1e6
    read_df_fold_change = np.log2(read_rpm_df.divide(read_rpm_df.input, axis=0)).drop(columns="input").sort_values("control").reset_index().rename(columns={"index": "chrom_coord"})
    # add starrpeaker annotations
    read_df_fold_change_with_annots = get_bed_annotations(read_df_fold_change, peak_files, peak_file_libname, "peak")
    # add diff_activity information
    diff_act_df = pd.concat([read_diff_output(df, dl) for df,dl in zip(diff_act_files, diff_act_libname)], axis=1)
    read_df_fold_change_with_annots = read_df_fold_change_with_annots.merge(diff_act_df, left_on="chrom_coord", right_index=True)
    return read_df_fold_change_with_annots


#################
# diff act file #
#################

def save_diff_act_file(df, store_dir, libn, activity_type):
    # save csv file
    save_file = get_diff_act_store_file(store_dir, libn, activity_type, "csv")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    df.to_csv(save_file, index=False)
    return


################
# multiprocess #
################

def run_singleargs_pool_job(pool_function, pool_iter):
    pool = mp.Pool(len(pool_iter))
    results = pool.map(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results

def run_multiargs_pool_job(pool_function, pool_iter):
    pool = mp.Pool(len(pool_iter))
    results = pool.starmap(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results
