import pandas as pd
import os

from cuppa.constants import MOCK_DATA_DIR
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.performance.performance_stats import PerformanceStats
from cuppa.visualization.visualization import CuppaVisData


class MockInputData:
    single_sample_tsv: str = os.path.join(MOCK_DATA_DIR, "input_data/COLO829v003T.cuppa_data.tsv.gz")
    cohort_dir: str = os.path.join(MOCK_DATA_DIR, "input_data/cohort/")


class MockTrainingData:

    path_features = os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz")
    path_metadata = os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv")
    path_fusion_overrides = os.path.join(MOCK_DATA_DIR, "training_data/fusion_overrides.tsv")

    X: pd.DataFrame = pd.read_table(path_features, index_col=0)
    metadata: pd.DataFrame = pd.read_table(path_metadata, index_col=0)
    y: pd.Series = metadata["CancerSubtype"]
    y_split: pd.Series = y + "__" + metadata["RnaReadLength"].astype(str)


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

    path_predictions: str = os.path.join(MOCK_DATA_DIR, "cv_output/predictions.tsv.gz")
    path_pred_summ: str = os.path.join(MOCK_DATA_DIR, "cv_output/pred_summ.tsv")
    path_performance: str = os.path.join(MOCK_DATA_DIR, "cv_output/performance.tsv")

    predictions = CuppaPrediction.from_tsv(path_predictions)
    pred_summ = CuppaPredSummary.from_tsv(path_pred_summ)
    performance = PerformanceStats.from_tsv(path_performance)
    probs_per_clf: dict[str, pd.DataFrame] = _get_probs_per_clf(predictions)


class MockCuppaClassifier:

    path_classifier: str = os.path.join(MOCK_DATA_DIR, "final_model/cuppa_classifier.pickle.gz")
    path_predictions: str = os.path.join(MOCK_DATA_DIR, "final_model/predictions.tsv.gz")

    cuppa_classifier: CuppaClassifier = CuppaClassifier.from_file(path_classifier)
    cuppa_classifier.cv_performance = MockCvOutput.performance

    predictions: CuppaPrediction = CuppaPrediction.from_tsv(path_predictions)
    probs_per_clf: dict[str, pd.DataFrame] = _get_probs_per_clf(predictions)


class MockVisData:

    path_vis_data: str = os.path.join(MOCK_DATA_DIR, "visualization/COLO829v003T.cuppa.vis_data.with_dummy_rna_preds.tsv")
    vis_data: CuppaVisData = CuppaVisData.from_tsv(path_vis_data)

