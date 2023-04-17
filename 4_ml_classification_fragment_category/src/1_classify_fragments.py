import os
import argparse
import numpy as np
import pandas as pd
import peaker as pk


def get_test_scores(cv_results):
    print(cv_results)
    best_estimator_idx = cv_results["rank_test_average_precision"][0] - 1 # index is 1-based
    aps = cv_results["mean_test_average_precision"][best_estimator_idx]
    prec = cv_results["mean_test_precision"][best_estimator_idx]
    rec = cv_results["mean_test_recall"][best_estimator_idx]
    acc = cv_results["mean_test_balanced_accuracy"][best_estimator_idx]
    f1 = cv_results["mean_test_f1"][best_estimator_idx]
    roc = cv_results["mean_test_roc_auc"][best_estimator_idx]
    return aps, prec, rec, acc, f1, roc

def create_peaker_compatible_input_df(file_dir, lib, file_base):
    file_path = os.path.join(file_dir, lib, f"{file_base}.bed")
    file_df = pd.read_csv(file_path, sep="\t", header=None, usecols=[3], names=["unique_id"])
    file_df["label"] = file_base
    return file_df.set_index("unique_id")

def create_peaker_compatible_inputs(file_dir, lib, file1_base, file2_base, store_dir):
    label_dict = {file1_base: 1, file2_base: 0}
    file1_df = create_peaker_compatible_input_df(file_dir, lib, file1_base)
    file2_df = create_peaker_compatible_input_df(file_dir, lib, file2_base)
    peaker_df = pd.concat([file1_df, file2_df], axis=0)
    store_path = os.path.join(store_dir, lib, f"{file1_base}_vs_{file2_base}", "labels.csv")
    os.makedirs(os.path.dirname(store_path), exist_ok=True)
    peaker_df.to_csv(store_path, index=False)
    return label_dict, peaker_df


def classify_enhancers(homer_pickle, class_file_dir, lib, class_type1, class_type2, save_dir, nthreads):
    # create peaker compatible inputs
    ld, df_labels = create_peaker_compatible_inputs(class_file_dir, lib, class_type1, class_type2, save_dir)
    # run peaker 
    bpk = pk.Biopeaker([df_labels], label_dict=ld, nthreads=int(nthreads))
    bpk.homer_featurize(homer_pickle=homer_pickle)
    cd,lm,gp,tsc = bpk.peaker()
    # get scores
    aps, prec, rec, acc, f1, roc = get_test_scores(gp.cv_results_)
    # save scores on validation
    save_test_score_file = os.path.join(save_dir, lib, f"{class_type1}_vs_{class_type2}", "scores.csv")
    os.makedirs(os.path.dirname(save_test_score_file), exist_ok=True)
    with open(save_test_score_file, "w") as f:
        f.write("APS,Precision,Recall,Accuracy,F1_score,ROCAUC\n")
        f.write(f"{aps},{prec},{rec},{acc},{f1},{roc}\n")
    # save contributing features
    save_feat_file = os.path.join(save_dir, lib, f"{class_type1}_vs_{class_type2}", "features.csv")
    os.makedirs(os.path.dirname(save_feat_file), exist_ok=True)
    ldr = {v:k for k,v in ld.items()}
    with open(save_feat_file, "w") as f:
        for tf, val in cd[ldr[1]]:
            f.write(f"{tf},{val}\n")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq enhancer classification")
    parser.add_argument("homer_pickle", type=str, help="the pickle file with features for the enhancers")
    parser.add_argument("lib", type=str, help="library shortform as given in the meta file")
    parser.add_argument("fragment_class_dir", type=str, help="Dir where the fragments classified by their activity is stored")
    parser.add_argument("class1", type=str, help="Fragment Class 1")
    parser.add_argument("class2", type=str, help="Fragment Class 2")
    parser.add_argument("store_dir", type=str, help="Output dir where results will be stored")
    parser.add_argument("--nthreads", type=int, help="Number of cores", default=8)

    cli_args = parser.parse_args()

    classify_enhancers(
        cli_args.homer_pickle,
        cli_args.fragment_class_dir,
        cli_args.lib,
        cli_args.class1,
        cli_args.class2,
        cli_args.store_dir,
        cli_args.nthreads
    )
