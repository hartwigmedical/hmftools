import pandas as pd
import os
import joblib

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.performance.performance_stats import PerformanceStats
from cuppa.misc.cached_class_property import cached_class_property
from cuppa.constants import MOCK_DATA_DIR

"""
@cached_class_property is used throughout this module so that data is only loaded when actually needed (lazy loading) 
"""

class MockInputData:
    path_tsv_new_format_colo = os.path.join(MOCK_DATA_DIR, "input_data/COLO829T.cuppa_data.tsv.gz")
    path_tsv_new_format_prostate = os.path.join(MOCK_DATA_DIR, "input_data/prostate_sample.cuppa_data.tsv.gz")


class MockTrainingData:

    path_features = os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz")
    path_metadata = os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv")

    @cached_class_property
    def X(self) -> pd.DataFrame:
        return pd.read_table(self.path_features, index_col=0)

    @cached_class_property
    def metadata(self) -> pd.DataFrame:
        return pd.read_table(self.path_metadata, index_col=0)

    @cached_class_property
    def y(self) -> pd.Series:
        return self.metadata["CancerSubtype"]

    @cached_class_property
    def y_split(self) -> pd.Series:
        return self.y + "__" + self.metadata["RnaReadLength"].astype(str)


class MockCvOutput:

    path_predictions = os.path.join(MOCK_DATA_DIR, "training_output/cv/predictions.tsv")
    path_pred_summ = os.path.join(MOCK_DATA_DIR, "training_output/cv/pred_summ.tsv")
    path_performance = os.path.join(MOCK_DATA_DIR, "training_output/cv/performance.tsv")

    path_predictions_for_vis = os.path.join(MOCK_DATA_DIR, "visualization/predictions.tsv.gz")

    @cached_class_property
    def predictions(self) -> CuppaPrediction:
        return CuppaPrediction.from_tsv(self.path_predictions)

    @cached_class_property
    def pred_summ(self) -> CuppaPredSummary:
        return CuppaPredSummary.from_tsv(self.path_pred_summ)

    @cached_class_property
    def performance(self) -> pd.DataFrame:
        return PerformanceStats.from_tsv(self.path_performance)

    @cached_class_property
    def predictions_for_vis(self) -> CuppaPrediction:
        return CuppaPrediction.from_tsv(self.path_predictions_for_vis)


class MockTrainingOutput:

    path_cuppa_classifier = os.path.join(MOCK_DATA_DIR, "training_output/cuppa_classifier.pickle.gz")

    @cached_class_property
    def cuppa_classifier(self) -> CuppaClassifier:
        return joblib.load(self.path_cuppa_classifier)


class MockProbsFromFitTransform:

    ## Sub-classifiers --------------------------------
    @cached_class_property
    def gen_pos(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/gen_pos.txt"), index_col=0)

    @cached_class_property
    def snv96(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/snv96.txt"), index_col=0)

    @cached_class_property
    def event(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/event.txt"), index_col=0)

    @cached_class_property
    def gene_exp(self) -> pd.DataFrame:
        probs = pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/gene_exp.txt"), index_col=0)
        probs = probs.dropna(axis=0, how="all").dropna(axis=1, how="all")
        return probs

    @cached_class_property
    def alt_sj(self) -> pd.DataFrame:
        probs = pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/alt_sj.txt"), index_col=0)
        probs = probs.dropna(axis=0, how="all").dropna(axis=1, how="all")
        return probs

    ## Meta-classifiers --------------------------------
    @cached_class_property
    def dna_combined(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/dna_combined.txt"), index_col=0)

    @cached_class_property
    def rna_combined(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/rna_combined.txt"), index_col=0)

    @cached_class_property
    def combined(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_fit_transform/combined.txt"), index_col=0)


class MockProbsPreCal:
    @cached_class_property
    def _df(self) -> pd.DataFrame:
        return pd.read_table(os.path.join(MOCK_DATA_DIR, "probs_precalibration/probs_precalibration.tsv"), index_col=0)

    @cached_class_property
    def y(self) -> pd.Series:
        return self._df["actual_class"]

    @cached_class_property
    def X(self) -> pd.DataFrame:
        return self._df.drop("actual_class", axis=1)



