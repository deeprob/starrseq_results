import os
import argparse
import pandas as pd
import numpy as np


def get_samples_file(samples_dir, lib, pred_type):
    return os.path.join(samples_dir, lib, pred_type, "samples.csv")

def main(
    samples_dir, libs, pred_type, valid_lib, test_lib, meta_exp_file, save_dir,
    rna_thresh=1
    ):
    rna_df = pd.read_csv(meta_exp_file, index_col=[0], usecols=["gene_name"]+libs)
    sel_df = rna_df.loc[(rna_df>rna_thresh).any(axis=1)]
    sel_df = sel_df.groupby("gene_name").mean()
    sel_df = (sel_df - sel_df.mean())/sel_df.std()

    print(pred_type)
    meta_samples_df = pd.DataFrame()
    meta_feat_df = pd.DataFrame()
    for lib in libs:
        print(lib)
        samples_file = get_samples_file(samples_dir, lib, pred_type)
        samples_df = pd.read_csv(samples_file, usecols=["chrm", "start", "end", "label"])
        feat_df = pd.DataFrame(np.repeat(sel_df[lib].values.reshape(1, -1), len(samples_df), axis=0), columns=list(sel_df[lib].index))
        feat_df.insert(0, "chrom_coord", samples_df.chrm + "_" + samples_df.start.astype(str) + "_" + samples_df.end.astype(str))
        if lib == test_lib:
            samples_df["split"] = "test"
        elif lib == valid_lib:
            samples_df["split"] = "valid"
        else:
            samples_df["split"] = "train"
        meta_samples_df = meta_samples_df.append(samples_df)
        meta_feat_df = meta_feat_df.append(feat_df)
    sample_save_dir = os.path.join(save_dir, pred_type)
    os.makedirs(sample_save_dir, exist_ok=True)
    sample_save_file = os.path.join(sample_save_dir, "samples.h5")
    feat_save_file = os.path.join(sample_save_dir, "rna_features.h5")
    meta_samples_df.to_hdf(sample_save_file, index=False, key="samples")
    meta_feat_df.to_hdf(feat_save_file, index=False, key="features")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq enhancer classification")
    parser.add_argument("--samples_dir", type=str, help="File where samples with their labels are stored")
    parser.add_argument("--libs", type=str, help="library names to be used to create samples", nargs="+")
    parser.add_argument("--pred_type", type=str, help="The type of prediction")
    parser.add_argument("--valid_lib", type=str, help="The validation library name")
    parser.add_argument("--test_lib", type=str, help="The test library name")
    parser.add_argument("--meta_exp_file", type=str, help="The meta expression filepath")
    parser.add_argument("--save_dir", type=str, help="Dir where the file will be stored")
    parser.add_argument("--rnaseq_thresh", type=float, help="Dir where the file will be stored", default=1.0)

    cli_args = parser.parse_args()

    main(
        cli_args.samples_dir, cli_args.libs, cli_args.pred_type, 
        cli_args.valid_lib, cli_args.test_lib, cli_args.meta_exp_file, 
        cli_args.save_dir, cli_args.rnaseq_thresh)
