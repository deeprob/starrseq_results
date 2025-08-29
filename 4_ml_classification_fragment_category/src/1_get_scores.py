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

def save_test_scores(pred_save_file, score_save_file):
    aps, prec, rec, acc, f1, roc = get_test_scores(pred_save_file)
    # save scores on test
    os.makedirs(os.path.dirname(score_save_file), exist_ok=True)
    with open(score_save_file, "w") as f:
        f.write("APS,Precision,Recall,Accuracy,F1_score,ROCAUC\n")
        f.write(f"{aps},{prec},{rec},{acc},{f1},{roc}\n")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq enhancer classification")
    parser.add_argument("pred_filepath", type=str, help="File where model predictions are stored")
    parser.add_argument("score_filepath", type=str, help="File where model scores will be stored")

    cli_args = parser.parse_args()

    save_test_scores(cli_args.pred_filepath, cli_args.score_filepath)