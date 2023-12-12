import pandas as pd
import os
import shutil
from sklearn.compose import make_column_selector

from cuppa.constants import MOCK_DATA_DIR
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction
from cuppa.performance.performance_stats import PerformanceStats
from cuppa.misc.stat_tests import chi2_test
from cuppa.runners import TrainingRunner


def make_training_data():
    input_dir = "/resources/model/training_data"

    ## Load data --------------------------------
    ## TODO: Refactor this section. DataLoader is deprecated and now deleted
    data_loader = DataLoader(data_dir=input_dir)

    data_loader.load_features_dna(excl_chroms="chrY", genome_version=37)
    data_loader.load_features_rna()
    data_loader.load_metadata()

    ## Sample selection --------------------------------
    sample_selection = pd.read_excel(os.path.join(MOCK_DATA_DIR,"sample_selection.xlsx"))
    sample_ids = sample_selection["sample_id"]
    sample_ids_anon = sample_selection["actual_class"].str[0:3] + "_" + sample_selection["number"].astype(str)

    metadata = data_loader.metadata.loc[sample_ids, ["CancerType", "CancerSubtype", "RnaReadLength"]]

    features = data_loader.get_features()
    features = features.loc[sample_ids]

    ## Anonymize sample ids --------------------------------
    features.index = sample_ids_anon
    metadata.index = sample_ids_anon

    ## Subset for top RNA features --------------------------------
    def select_k_best_features(X, y, k=1000) -> pd.Series:
        #X=features[make_column_selector("gene_exp")]
        #y=labels["CancerSubtype"]
        #k=1000

        X = X.dropna()
        y = y[X.index]

        chi2_out = chi2_test(X,y)
        sel_features = chi2_out["feature"][0:k].sort_values()

        return sel_features

    features_gene_exp = features[ select_k_best_features(features[make_column_selector("gene_exp")], metadata["CancerSubtype"]) ]
    features_alt_sj = features[ select_k_best_features(features[make_column_selector("alt_sj")], metadata["CancerSubtype"]) ]

    features_rna = pd.concat([features_gene_exp, features_alt_sj], axis=1)
    features_dna = features.loc[:,~features.columns.str.match("^(gene_exp|alt_sj)")]

    features = pd.concat([features_dna, features_rna], axis=1)

    ## Export --------------------------------
    output_dir = os.path.join(MOCK_DATA_DIR,"training_data")

    features.to_csv(os.path.join(output_dir,"features.tsv.gz"), sep="\t")
    metadata.to_csv(os.path.join(output_dir, "metadata.tsv"), sep="\t")

def make_training_output():
    ## Load input data --------------------------------
    X = pd.read_csv(os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz"), sep='\t', index_col=0)

    labels = pd.read_csv(os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv"), sep='\t', index_col=0)
    y = labels["CancerSubtype"]
    y_split = labels["CancerSubtype"] + "__" + labels["RnaReadLength"].astype(str)

    ## Train --------------------------------
    runner = TrainingRunner(
        #fusion_overrides_path=None,  ## Ignore fusion override
        input_dir="PLACEHOLDER",
        output_dir="PLACEHOLDER",
        n_jobs=1
    )

    ## Input data
    runner.X = X
    runner.y = y
    runner.y_split = y_split

    ## Paths
    output_dir = os.path.join(MOCK_DATA_DIR, "model")
    os.makedirs(output_dir, exist_ok=True)
    runner.output_dir = output_dir

    ## Train
    runner.cv_fit()
    runner.train_final_model()

    ## Write output --------------------------------
    cv_out_dir = os.path.join(MOCK_DATA_DIR, "training_output/cv")
    os.makedirs(cv_out_dir, exist_ok=True)

    runner.cv_predictions.to_tsv(os.path.join(cv_out_dir, "predictions.tsv"))
    runner.cv_pred_summ.to_tsv(os.path.join(cv_out_dir, "pred_summ.tsv"))
    runner.cv_performance.to_csv(os.path.join(cv_out_dir, "performance.tsv"), sep="\t", index=False)

    runner.cuppa_classifier.to_file(os.path.join(MOCK_DATA_DIR, "training_output/cuppa_classifier.pickle.gz"))

    ## Delete training output folder
    shutil.rmtree(output_dir)

def make_probs_intermediate():
    X = pd.read_table(os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz"), index_col=0)
    labels = pd.read_table(os.path.join(MOCK_DATA_DIR, "training_data/labels.tsv"), index_col=0)
    y = labels["CancerSubtype"]

    cuppa_classifier = CuppaClassifier()
    cuppa_classifier.fit(X, y)

    probs = cuppa_classifier.predict_proba(X)

    out_dir = os.path.join(MOCK_DATA_DIR, "probs_fit_transform/")
    clf_names = probs.index.get_level_values("clf_name").unique()
    for clf_name in clf_names:
        probs_clf = probs[ probs.index.get_level_values("clf_name")==clf_name ].copy()
        probs_clf.index = probs_clf.index.get_level_values("sample_id")
        probs_clf.to_csv(os.path.join(out_dir, clf_name+".txt"), sep="\t")


def get_predictions_for_vis():

    ## Load data
    cv_report_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/04-NET_as_prefix/cv/report/"
    predictions = CuppaPrediction.from_tsv(os.path.join(cv_report_dir, "predictions.tsv.gz"))
    cv_performance = PerformanceStats.from_tsv(os.path.join(cv_report_dir, "perf.tsv"))

    ## Select samples
    sample_selection = pd.read_excel(os.path.join(MOCK_DATA_DIR, "sample_selection.xlsx"), sheet_name="sample_selection.vis")
    predictions_subset = predictions.get_samples(sample_selection["sample_id"])

    ## Anonymize sample_ids
    index = predictions_subset.index.to_frame(index=False)
    sample_id_mappings = pd.Series(sample_selection["sample_id_anon"].values, sample_selection["sample_id"])
    index["sample_id"] = sample_id_mappings[index["sample_id"]].values
    predictions_subset.index = pd.MultiIndex.from_frame(index)

    ## Add cv performance
    predictions_subset = predictions_subset.add_cv_performance(cv_performance)

    predictions_subset.to_tsv(os.path.join(MOCK_DATA_DIR, "visualization/predictions.tsv.gz"))
