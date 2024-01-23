import os
import shutil
import sys
import tempfile

import pandas as pd

from cuppa.tests.mock_data import MockTrainingData, MockCuppaClassifier, MockInputData, MockCvOutput
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.constants import DEFAULT_CUPPA_CLASSIFIER_PATH
from cuppa.runners import PredictionRunner, TrainingRunner, RunnerArgParser
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary


class TestRunnerArgParser:
    def test_can_create_prediction_runner_from_required_command_line_args(self):

        sys.argv = [
            "example_module.py",
            "--features_path=/path/to/features/",
            "--output_dir=/path/to/output/",
            "--classifier_path=/path/to/cuppa_classifier.pickle.gz"
        ]

        arg_parser = RunnerArgParser()
        kwargs = arg_parser.get_kwargs_predict()

        runner = PredictionRunner(**kwargs)
        assert True

    def test_can_create_training_runner_from_required_command_line_args(self):
        sys.argv = [
            "example_module.py",
            "--features_path=/path/to/features/",
            "--metadata_path=/path/to/metadata_file",
            "--output_dir=/path/to/output/",
        ]

        arg_parser = RunnerArgParser()
        kwargs = arg_parser.get_kwargs_train()

        runner = TrainingRunner(**kwargs)
        assert True


class TestPredictionRunner:

    def test_predict_with_real_classifier_and_sample_gives_correct_results_and_writes_output(self):

        output_dir = os.path.join(tempfile.gettempdir(), "pycuppa_prediction_run_test")
        os.makedirs(output_dir, exist_ok=True)

        runner = PredictionRunner(
            features_path=MockInputData.path_tsv_new_format_prostate,
            sample_id="TEST_SAMPLE",
            classifier_path=DEFAULT_CUPPA_CLASSIFIER_PATH,
            output_dir=output_dir,
            using_old_features_format=False
        )

        runner.run()

        combined_prediction = runner.pred_summ.query("clf_name=='combined'")
        assert combined_prediction["pred_prob_1"].round(4).iloc[0] == 0.9968
        assert combined_prediction["pred_class_1"].iloc[0] == "Prostate"

        rna_combined_prediction = runner.pred_summ.query("clf_name=='rna_combined'")
        assert rna_combined_prediction["pred_prob_1"].round(4).iloc[0] == 0.9953
        assert rna_combined_prediction["pred_class_1"].iloc[0] == "Prostate"

        assert os.path.exists(runner.plot_path)
        assert os.path.exists(runner.vis_data_path)
        assert os.path.exists(runner.pred_summ_path)

        shutil.rmtree(output_dir)

    def test_fusion_override_alters_probabilities(self):

        output_dir = os.path.join(tempfile.gettempdir(), "pycuppa_prediction_run_test")
        os.makedirs(output_dir, exist_ok=True)

        runner = PredictionRunner(
            features_path=MockInputData.path_tsv_new_format_prostate,
            sample_id="TEST_SAMPLE",
            classifier_path=DEFAULT_CUPPA_CLASSIFIER_PATH,
            output_dir=output_dir,
            using_old_features_format=False
        )

        runner.run()

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


class TestTrainingRunner:

    # from cuppa.tests.test_runners import TestTrainingRunner
    # self = TestTrainingRunner

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

    def _run_individual_steps_on_real_data(self):

        input_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/"

        runner = TrainingRunner(
            features_path=input_dir,
            output_dir="/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/07-pip_env_2",
            metadata_path=os.path.join(input_dir, "cup_ref_sample_data.csv"),
            log_to_file=True,
            using_old_features_format=True,
            n_jobs=5
        )

        runner.load_data()
        runner.cv_fit()
        runner.make_cv_report()

        runner.train_final_model()
        runner.export_final_model()
