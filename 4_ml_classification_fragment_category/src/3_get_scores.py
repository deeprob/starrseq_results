import os
import argparse
import pandas as pd
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score, f1_score, balanced_accuracy_score, accuracy_score, precision_score, recall_score


def get_test_scores(pred_save_file, thresh=0.5):
    df = pd.read_csv(pred_save_file, header=None, names=["ypred", "ytarget", "chrm", "start", "end"])
    aps = average_precision_score(df.ytarget, df.ypred)
    roc = roc_auc_score(df.ytarget, df.ypred)
    prec = precision_score(df.ytarget, (df.ypred>thresh).astype(int))
    rec = recall_score(df.ytarget, (df.ypred>thresh).astype(int))
    acc = balanced_accuracy_score(df.ytarget, (df.ypred>thresh).astype(int))
    f1 = f1_score(df.ytarget, (df.ypred>thresh).astype(int))
    return aps, prec, rec, acc, f1, roc

def save_test_scores(proj_dir, lib, model, class1, class2):
    save_file = os.path.join(proj_dir, lib, f"{class1}_vs_{class2}", model, f"{model}.csv.gz")
    aps, prec, rec, acc, f1, roc = get_test_scores(save_file)
    # save scores on test
    save_test_score_file = os.path.join(proj_dir, lib, f"{class1}_vs_{class2}", model, "scores.csv")
    os.makedirs(os.path.dirname(save_test_score_file), exist_ok=True)
    with open(save_test_score_file, "w") as f:
        f.write("APS,Precision,Recall,Accuracy,F1_score,ROCAUC\n")
        f.write(f"{aps},{prec},{rec},{acc},{f1},{roc}\n")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq enhancer classification")
    parser.add_argument("lib", type=str, help="library shortform")
    parser.add_argument("model", type=str, help="the model used to classify")
    parser.add_argument("proj_dir", type=str, help="project dir where the fragments classified by their activity is stored")
    parser.add_argument("class1", type=str, help="Fragment Class 1")
    parser.add_argument("class2", type=str, help="Fragment Class 2")

    cli_args = parser.parse_args()

    save_test_scores(cli_args.proj_dir, cli_args.lib, cli_args.model, cli_args.class1, cli_args.class2)
