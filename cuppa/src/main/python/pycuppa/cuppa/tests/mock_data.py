import pandas as pd
import os

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.performance.performance_stats import PerformanceStats
from cuppa.misc.cached_class_property import cached_class_property
from cuppa.constants import MOCK_DATA_DIR
from cuppa.visualization.visualization import CuppaVisData

"""
@cached_class_property is used throughout this module so that data is: 
- only loaded when actually needed (lazy loading)
- immutable 
"""


class MockInputData:
    path_tsv_new_format = os.path.join(MOCK_DATA_DIR, "input_data/new_format/COLO829v003T.cuppa_data.tsv.gz")
    dir_old_format = os.path.join(MOCK_DATA_DIR, "input_data/old_format")


class MockTrainingData:

    path_features = os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz")
    path_metadata = os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv")
    path_fusion_overrides = os.path.join(MOCK_DATA_DIR, "training_data/fusion_overrides.tsv")

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


def _get_probs_per_clf(predictions: CuppaPrediction) -> dict[str, pd.DataFrame]:

    probs = predictions[predictions.index.get_level_values("data_type") == "prob"]
    sample_order = MockTrainingData.y.index
    probs = probs.loc[sample_order]

    clf_names = probs.index.get_level_values("clf_name").unique()
    probs_per_clf = {}
    for clf_name in clf_names:
        probs_i = probs[probs.index.get_level_values("clf_name") == clf_name].copy()
        probs_i.index = probs_i.index.get_level_values("sample_id")
        probs_i.columns.name = None
        probs_per_clf[clf_name] = probs_i

    return probs_per_clf


class MockCvOutput:

    path_predictions = os.path.join(MOCK_DATA_DIR, "cv_output/predictions.tsv.gz")
    path_pred_summ = os.path.join(MOCK_DATA_DIR, "cv_output/pred_summ.tsv")
    path_performance = os.path.join(MOCK_DATA_DIR, "cv_output/performance.tsv")

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
    def probs_per_clf(self) -> dict[str, pd.DataFrame]:
        return _get_probs_per_clf(self.predictions)


class MockCuppaClassifier:

    path_classifier = os.path.join(MOCK_DATA_DIR, "final_model/cuppa_classifier.pickle.gz")
    path_predictions = os.path.join(MOCK_DATA_DIR, "final_model/predictions.tsv.gz")

    @cached_class_property
    def cuppa_classifier(self) -> CuppaClassifier:
        cuppa_classifier = CuppaClassifier.from_file(self.path_classifier)
        cuppa_classifier.cv_performance = MockCvOutput.performance
        return cuppa_classifier

    @cached_class_property
    def predictions(self) -> CuppaPrediction:
        return CuppaPrediction.from_tsv(self.path_predictions)

    @cached_class_property
    def probs_per_clf(self) -> dict[str, pd.DataFrame]:
        return _get_probs_per_clf(self.predictions)


class MockVisData:

    path_vis_data = os.path.join(MOCK_DATA_DIR, "visualization/COLO829v003T.cuppa.vis_data.with_dummy_rna_preds.tsv")

    @cached_class_property
    def vis_data(self) -> CuppaVisData:
        return CuppaVisData.from_tsv(self.path_vis_data)
