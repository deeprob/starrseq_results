import os
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, average_precision_score, f1_score, roc_auc_score

class HomerVectorizer:
    """
    Vectorizes chromosomal locations by scanning them using Homer and its motif database
    """
    def __init__(self, roi_bed, genome_fasta, homer_pwm_motifs, homer_outdir, threads=32):
        self.genome = genome_fasta
        self.homer_pwms = homer_pwm_motifs
        self.roi = roi_bed
        self.homer_outdir = homer_outdir
        self.threads = threads

        self.homer_outfile = self._get_homer_outfile()
        self.homer_roi = self._get_roi_homer()
        self.logfile = self._get_logfile()
        pass

    def _get_homer_outfile(self):
        return os.path.join(self.homer_outdir, "motif_odds.tsv")

    def _get_roi_homer(self):
        return os.path.join(self.homer_outdir, "tmp_roi_homer.bed")

    def _get_logfile(self):
        return os.path.join(self.homer_outdir, "homer_scan.log")

    def _process_homer(self, df_row):
        chrname = df_row.chrm
        start = df_row.start
        end = df_row.end
        peak_name = f"{chrname}_{start}_{end}"
        irrelevant_col = 0
        strand = "."
        return pd.Series({"3":peak_name, "4":irrelevant_col, "5": strand}) 

    def _create_homer_compatible_roi(self):
        """
        Converts an roi file with chromosomal coordinates to a homer compatible one
        """
        df_roi = pd.read_csv(self.roi, usecols=[0,1,2], sep="\t", header=None).drop_duplicates()
        df_roi.rename(columns={0: "chrm", 1: "start", 2: "end"}, inplace=True)
        df_roi.loc[:, ["chrm", "start", "end"]].merge(df_roi.apply(self._process_homer, axis=1), left_index=True, right_index=True).to_csv(self.homer_roi, index=False, header=None, sep="\t")
        return

    def _pwm_scan_homer(self):

        cmd = [
            "findMotifsGenome.pl", self.homer_roi, self.genome, self.homer_outdir, 
            "-find", self.homer_pwms, "-p", str(self.threads), "-size", "given", 
            ]
        with open(self.logfile, "w") as lf:
            with open(self.homer_outfile, "w") as of:
                results = subprocess.run(cmd, stdout=of, stderr=lf)
        return results

    def _parse_homer_outfile(self):
        motif_df = pd.read_csv(self.homer_outfile, sep="\t")
        # pivot the table for easier functioning, 
        # also for all motifs present more than once in a region scores will be summed 
        motif_df = motif_df.pivot_table('MotifScore', ["PositionID"], ["Motif Name", "Strand"], aggfunc=np.sum, fill_value=0.0)
        # fix column names and their levels after pivot table
        motif_df.columns = [f'{i}|{j}' if j != '' else f'{i}' for i,j in motif_df.columns]
        return motif_df

    def featurize(self):
        # making sure there is a directory on path
        os.makedirs(self.homer_outdir, exist_ok=True)
        # creating homer compatible tmp roi file
        self._create_homer_compatible_roi()
        # running homer
        self._pwm_scan_homer()
        # parsing the homer outfile
        self.motif_df = self._parse_homer_outfile()
        return

    def store_features_to_pickle(self, pickle_path):
        self.motif_df.to_pickle(pickle_path)
        return

    @classmethod
    def load_features_from_pickle(self, pickle_path):
        return pd.read_pickle(pickle_path)


class Biopeaker:

    def __init__(self, labelized_dfs, label_dict=dict(), nthreads=-1):
        """
        labelized_df: A dataframe or pandas series object with chromosomal coordinates as index and their labels as the label column
        """
        self.featurized_df = pd.DataFrame()
        self.labelized_dfs = labelized_dfs
        if not label_dict:
            uniq_vals = labelized_dfs[0].label.unique()
            label_dict = dict(zip(uniq_vals, range(len(uniq_vals))))
        self.label_dict = label_dict
        self.nthreads = nthreads
        pass

    ##############
    # featurizer #
    ##############

    def homer_featurize(self, **kwargs):

        if "homer_pickle" in kwargs:
            homer_pickle = kwargs["homer_pickle"]
            if os.path.exists(homer_pickle):
                self.featurized_df = HomerVectorizer.load_features_from_pickle(homer_pickle)
                print("Loaded pickle ... ")
            else:
                raise IOError("FileNotFound: Pickle file does not exist!")
        else:
            # assert that the required arguments for homer featurizer is present
            roi_bed = kwargs["roi_bed"]
            genome_fasta = kwargs["genome_fasta"]
            homer_pwm_motifs = kwargs["homer_pwm_motifs"]
            homer_outdir = kwargs["homer_outdir"]
            # run homer to vectorize regions
            hv = HomerVectorizer(roi_bed, genome_fasta, homer_pwm_motifs, homer_outdir)
            hv.featurize()
            homer_pickle = os.path.join(homer_outdir, "motif_features.pkl")
            hv.store_features_to_pickle(homer_pickle)
            print(f"Homer pickle stored in {homer_pickle}, Can be loaded directly from path next time!!!")
            self.featurized_df = HomerVectorizer.load_features_from_pickle(homer_pickle)

        return

    #############
    # labelizer #
    #############

    def _create_numeric_labels(self, df):
        df["label"] = df["label"].map(self.label_dict)
        return df


    ############################
    # preprocessing model data #
    ############################

    def _create_processed_df(self, labelized_df):
        """
        A function that creates a final processed dataframe which
        has chromosomal locations as index, features as the first N-1 columns 
        and labels as the final columns 
        """
        labelized_df = self._create_numeric_labels(labelized_df)
        motif_df_processed = self.featurized_df.merge(labelized_df, left_index=True,  right_index=True)
        return motif_df_processed

    def _get_processed_data(self, labelized_df):
        motif_df_processed = self._create_processed_df(labelized_df)
        X = motif_df_processed.iloc[:, :-1].values
        y = motif_df_processed.iloc[:, -1].values
        return X, y

    def _get_train_test(self):
        # check if labellized
        if not self.labelized_dfs:
            raise ValueError("Regions are not labelized: Please labelize before calling peaker")

        X_train, y_train = self._get_processed_data(self.labelized_dfs[0])

        test_data = []
        if len(self.labelized_dfs)>1:
            test_data = [self._get_processed_data(l_df) for l_df in self.labelized_dfs[1:]]    
        return X_train, y_train, test_data

    ##########
    # models #
    ##########

    def _get_linear_model(self, X_train, y_train):
        cw = {c:1/freq for c,freq in  dict(zip(*np.unique(y_train, return_counts=True))).items()}
        linear_baseline = make_pipeline(
            StandardScaler(), 
            LogisticRegression(max_iter=1000, class_weight=cw, penalty="l2", solver="lbfgs")
            )

        parameters = {
            "logisticregression__C": [0.1, 1, 10, 100, 1000],
            "logisticregression__penalty": ["l2"]
        }
        grid_pipeline = GridSearchCV(
            linear_baseline, 
            parameters, 
            scoring=["average_precision", "precision", "recall", "balanced_accuracy", "f1", "roc_auc"], 
            n_jobs=self.nthreads,
            cv=5,
            refit="average_precision"
            )
        grid_pipeline.fit(X_train, y_train)
        return grid_pipeline.best_estimator_, grid_pipeline

    ############
    # evaluate #
    ############

    def _get_sklearn_model_preds(self, model, X_test):
        y_pred = model.predict(X_test)
        return y_pred

    def _get_sklearn_model_prob(self, model, X_test):
        y_prob = model.predict_proba(X_test)
        return y_prob[:, 1]

    def _get_evaluation_scores(self, y_test, y_pred, y_pred_prob):
        acc = accuracy_score(y_test, y_pred)
        bacc = balanced_accuracy_score(y_test, y_pred)
        pr = precision_score(y_test, y_pred)
        re = recall_score(y_test, y_pred)
        aps = average_precision_score(y_test, y_pred_prob)
        f1 = f1_score(y_test, y_pred)
        roc = roc_auc_score(y_test, y_pred_prob)
        return acc, bacc, pr, re, aps, f1, roc

    def _get_coeff_dict(self, linear_model):
        coeff_dict = {}
        label_dict_rev = {v: k for k,v in self.label_dict.items()}
        if len(self.label_dict.keys()) == 2:
            # binary classification
            # coefficient array of dimension (1, n_features)
            map_tuple = tuple(zip(list(self.featurized_df.columns),linear_model['logisticregression'].coef_[0]))
            # the coeffs correspond to the label assigned numerical 1 
            coeff_dict[label_dict_rev[1]] = sorted(map_tuple, key=lambda x:abs(x[1]), reverse=True)
        else:
            for num_val, label in self.label_dict.items():
                map_tuple = tuple(zip(list(self.featurized_df.columns),linear_model['logisticregression'].coef_[num_val, :]))
                coeff_dict[label] = sorted(map_tuple, key=lambda x:abs(x[1]), reverse=True)
        return coeff_dict

    ##########
    # peaker #
    ##########

    def peaker(self):

        # check if featurized
        if self.featurized_df.empty:
            raise ValueError("Regions are not featurized: Please featurize before calling peaker")
        
        # get train-valid-test data
        X_train, y_train, test_data = self._get_train_test()

        # train linear model
        linear_model, grid_pipeline = self._get_linear_model(X_train, y_train)

        # get the feature coeffs
        coeff_dict = self._get_coeff_dict(linear_model)

        # evaluate on test using average precision score
        test_scores = []
        for X_test, y_test in test_data:
            y_pred = self._get_sklearn_model_preds(linear_model, X_test)
            y_prob = self._get_sklearn_model_prob(linear_model, X_test)
            scores = self._get_evaluation_scores(y_test, y_pred, y_prob)
            test_scores.append(scores)

        return coeff_dict, linear_model, grid_pipeline, test_scores
