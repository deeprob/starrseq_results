import os
import subprocess
import pandas as pd
import pybedtools
import json
from argparse import Namespace
import itertools
import functools
import multiprocessing as mp
from pybedtools.featurefuncs import greater_than

# TODO: take tmp dir from arguments
pybedtools.helpers.set_tempdir("/data5/deepro/tmp/")

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
        genome_sizes = meta_dict["genome"]["chrom_sizes"],
        roi_file = meta_dict["roi"]["sorted"],
    )

    return args


####################
# Filename parsing #
####################

def get_depth_bed_filepaths(store_dir, lib_short, lib_prefix, lib_rep):
    return os.path.join(store_dir, lib_short, f"{lib_prefix}_{lib_rep}.bed")

def get_lib_depth_beds_filepaths(store_dir, lib_short, lib_prefix, lib_reps):
    depth_beds = [get_depth_bed_filepaths(store_dir, lib_short, lib_prefix, lib_rep) for lib_rep in lib_reps.split()]
    return depth_beds


##################
# create windows #
##################

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

def make_windows_faster(in_bed, out_bed, window_size=500, window_stride=50):
    """
    Break the ROIs into fragments of an user defined window size and stride
    """
    window = pybedtools.BedTool().window_maker(b=in_bed, w=window_size, s=window_stride)
    # only keep windows of length greater than w-s-1
    window_new_filtered = window.filter(greater_than, window_size - window_stride + 1)
    window_new_filtered.saveas(out_bed)
    return

#############
# ROI depth #
#############

def get_roi_depth(filtered_bam, roi_sorted_bed, bed_out):
    # bam = pybedtools.BedTool(filtered_bam)
    roi = pybedtools.BedTool(roi_sorted_bed)
    try:
        c = roi.coverage(filtered_bam, sorted=True)
    except Exception as e:
        c = roi.coverage(filtered_bam)
    os.makedirs(os.path.dirname(bed_out), exist_ok=True)
    c.moveto(bed_out)
    return

#########################
# diff activity helpers #
#########################

def get_rep_depth_files(depth_dir, lib_short, lib_prefix, lib_reps):
    depth_files = [os.path.join(depth_dir, lib_short, f"{'_'.join([lib_prefix, rep])}.bed") for rep in lib_reps.split()]
    return depth_files


def get_lib_diff_activity_peak_filepath(store_dir, lib_short, diff_activity_type):
    peak_filepath = os.path.join(
        store_dir, lib_short, f"{diff_activity_type}.bed"
        )
    return peak_filepath

def get_replicate_depth_ser(depth_file):
    depth_df =  pd.read_csv(depth_file, header=None, sep="\t")
    depth_df = depth_df.set_index([0,1,2])
    depth_df.index = depth_df.index.rename(["chrom", "start", "end"])
    depth_reads = depth_df.iloc[:, -4]
    return depth_reads

def get_replicate_wise_depth_df(depth_files, lib_name):
    rep_depth_sers = [get_replicate_depth_ser(df) for df in depth_files] # multi_args_pool_job(get_replicate_norm_depth, pool_iter)
    rep_depth_dfs = pd.concat(rep_depth_sers, axis=1)
    if lib_name[0].isdigit(): # check to see if libname starts with a digit, else deseq2 will raise an error because R is horrible
        rep_depth_dfs.columns = [f"X{lib_name}_R{i}" for i in range(1, rep_depth_dfs.shape[1] + 1)]    
    else:
        rep_depth_dfs.columns = [f"{lib_name}_R{i}" for i in range(1, rep_depth_dfs.shape[1] + 1)]
    return rep_depth_dfs

def save_deseq_in_df(depth_df, store_dir, lib_short):
    save_file = os.path.join(store_dir, lib_short, "deseq_in.csv")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    # convert the chromosomal coordinates to proper format
    depth_df = depth_df.reset_index().drop_duplicates(keep="first")
    depth_df["unique_id"] = depth_df.chrom + "_" + depth_df.start.astype("str") + "_" + depth_df.end.astype("str")
    depth_df = depth_df.drop(columns= ["chrom", "start", "end"]).set_index("unique_id")
    depth_df.to_csv(save_file, index=True)    
    return save_file

def save_design_matrix(colnames, conditions, savedir, lib_short):
    savefile = os.path.join(savedir, lib_short, "deseq_design.csv")
    with open(savefile, "w") as f:
        f.write(f",condition\n")
        for col,line in zip(colnames, conditions):
            f.write(",".join([col,line]))
            f.write("\n")
    return savefile

def get_deseq_out_file(store_dir, lib_short):
    save_file = os.path.join(store_dir, lib_short, "deseq_out.csv")
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    return save_file

def da_with_deseq(table_in, design_mat_in, table_out):
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/scripts/run_DESeq2.sh", 
        table_in, design_mat_in, table_out
        ]
    subprocess.run(cmd)
    return

def save_deseq_diff_activity_file(deseq_outfile, store_dir, lib_short):
    df = pd.read_csv(deseq_outfile).reset_index().rename(columns={"index": "unique_id"})
    df = pd.concat((df, df.unique_id.str.split("_", expand=True).rename(columns={0: "chrom", 1: "start", 2: "end"})), axis=1)
    df["strand"] = "."
    df = df.loc[:, ["chrom", "start", "end", "unique_id", "stat", "strand", "log2FoldChange", "baseMean", "lfcSE", "pvalue", "padj"]]

    df_induced = df.loc[((df.log2FoldChange>0) & (df.padj<0.01))]
    df_repressed = df.loc[((df.log2FoldChange<0) & (df.padj<0.01))]
    df_nonresponsive = df.loc[~(((df.log2FoldChange>0) & (df.padj<0.01)) | ((df.log2FoldChange<0) & (df.padj<0.01)))]

    assert len(df) == sum(list(map(len, [df_induced, df_repressed, df_nonresponsive])))

    induced_file = os.path.join(store_dir, lib_short, "induced.bed")
    repressed_file = os.path.join(store_dir, lib_short, "repressed.bed")
    nonresponsive_file = os.path.join(store_dir, lib_short, "nonresponsive.bed")

    df_induced.to_csv(induced_file, sep="\t", header=False, index=False)
    df_repressed.to_csv(repressed_file, sep="\t", header=False, index=False)
    df_nonresponsive.to_csv(nonresponsive_file, sep="\t", header=False, index=False)
    return

def get_diff_active_fragments_deseq(
    lib_depth_files,
    cc_depth_files, 
    lib_short,
    cc_short, 
    store_dir
    ):
    # get the depth file for each replicate
    cc_depth_df = get_replicate_wise_depth_df(cc_depth_files, cc_short)
    lib_depth_df = get_replicate_wise_depth_df(lib_depth_files, lib_short)
    combined_depth_df = pd.concat([cc_depth_df, lib_depth_df], axis=1)
    # save the depth df 
    deseq_in_file = save_deseq_in_df(combined_depth_df, store_dir, lib_short)
    # get conditions from the depth dfs
    conditions = [cc_short for i in range(cc_depth_df.shape[1])] + [lib_short for i in range(lib_depth_df.shape[1])] 
    # generate and save the design matrix
    design_mat_file = save_design_matrix(combined_depth_df.columns, conditions, store_dir, lib_short)
    # get deseq out table
    deseq_out = get_deseq_out_file(store_dir, lib_short)
    # run deseq
    da_with_deseq(deseq_in_file, design_mat_file, deseq_out)
    # save induced, repressed and non responsive regions
    save_deseq_diff_activity_file(deseq_out, store_dir, lib_short)
    return

################
# multiprocess #
################

def run_singleargs_pool_job(pool_function, pool_iter, threads=None):
    if not threads:
        threads = len(pool_iter)    
    pool = mp.Pool(threads)
    results = pool.map(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results

def run_multiargs_pool_job(pool_function, pool_iter, threads=None):
    if not threads:
        threads = len(pool_iter)
    pool = mp.Pool(threads)
    results = pool.starmap(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results
