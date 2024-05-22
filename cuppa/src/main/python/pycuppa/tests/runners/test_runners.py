import os
import tempfile
import pandas as pd

from tests.mock_data import MockTrainingData, MockCuppaClassifier, MockCvOutput
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.runners.prediction_runner import PredictionRunner
from cuppa.runners.training_runner import TrainingRunner
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary


class TestPredictionRunner:

    def test_predict_with_mock_classifier_and_data_gives_correct_results(self):

        runner = PredictionRunner(
            features_path="/PLACEHOLDER/PATH/",
            output_dir="/PLACEHOLDER/PATH/",
            classifier_path=MockCuppaClassifier.path_classifier,
        )

        ## Mock features is already in 'samples x features' matrix form and doesn't need to be parsed
        ## Therefore manually assign features into runner to skip the parsing steps
        runner.X = MockTrainingData.X

        runner.get_predictions()
        runner.get_pred_summ()

        assert isinstance(runner.predictions, CuppaPrediction)
        assert isinstance(runner.pred_summ, CuppaPredSummary)

        sample_prediction = runner.pred_summ.query("sample_id=='35_Prostate' & clf_name=='combined'")
        assert sample_prediction["pred_prob_1"].round(3).iloc[0] == 0.944
        assert sample_prediction["pred_class_1"].iloc[0] == "Prostate"

        sample_prediction = runner.pred_summ.query("sample_id=='66_Lung' & clf_name=='rna_combined'")
        assert sample_prediction["pred_prob_1"].round(3).iloc[0] == 0.857
        assert sample_prediction["pred_class_1"].iloc[0] == "Lung"

    def test_can_look_up_cross_validation_predictions(self):

        runner = PredictionRunner(
            features_path="/PLACEHOLDER/PATH/",
            output_dir="/PLACEHOLDER/PATH/",
            classifier_path=MockCuppaClassifier.path_classifier,
        )
        runner.X = MockTrainingData.X

        cv_samples = pd.Series(["0_Breast", "1_Breast"])
        cv_predictions = MockCvOutput.predictions.loc[cv_samples]

        non_cv_samples = runner.X.index[~runner.X.index.isin(cv_samples)]

        ## Predict without providing CV predictions
        runner.get_predictions()
        predictions_all = runner.predictions.copy()

        ## Predict, providing CV predictions
        runner.cv_predictions = cv_predictions
        runner.get_predictions()
        predictions_with_lookup = runner.predictions.copy()

        ## Ensure rows are aligned
        predictions_with_lookup = predictions_with_lookup.reindex(predictions_all.index)

        assert \
            predictions_with_lookup.loc[non_cv_samples]\
            .equals(predictions_all.loc[non_cv_samples])

        assert \
            predictions_with_lookup.loc[cv_samples]\
            .equals(cv_predictions.loc[cv_samples])\
            .__invert__()


class TestTrainingRunner:

    # from tests.runners.test_runners import TestTrainingRunner; self = TestTrainingRunner

    output_dir = os.path.join(tempfile.gettempdir(), "pycuppa_prediction_run_test")
    os.makedirs(output_dir, exist_ok=True)

    runner = TrainingRunner(
        features_path="/PLACEHOLDER/PATH/",
        metadata_path="/PLACEHOLDER/PATH/",
        output_dir=output_dir,

        fusion_overrides_path=None,
        n_jobs=1,
        cache_training=False,
    )

    runner.X = MockTrainingData.X
    runner.y = MockTrainingData.y
    runner.y_split = MockTrainingData.y_split

    def test_cv_on_mock_data_gives_same_performance(self):
        self.runner.cv_fit()
        self.runner.get_cv_predictions()
        self.runner.get_cv_pred_summ()
        self.runner.get_cv_performance()

        actual_performance = self.runner.cv_performance.copy()
        expected_performance = MockCvOutput.performance.copy()

        ## Align rows
        actual_performance = actual_performance.set_index(["class", "clf_name"])
        expected_performance = expected_performance.set_index(["class", "clf_name"])

        comparison = pd.concat(
            dict(actual=actual_performance, expected=expected_performance),
            axis=1
        )

        ## Fill NAs with zero so that equal comparison works
        comparison = comparison.fillna(0).round(3)

        assert all(comparison[("actual","recall")] == comparison[("expected","recall")])
        assert all(comparison[("actual","precision")] == comparison[("expected","precision")])

    def test_train_final_model_successful(self):
        self.runner.train_final_model()
        assert isinstance(self.runner.cuppa_classifier, CuppaClassifier)
