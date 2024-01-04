import os
import shutil
import sys
import tempfile

import pandas as pd

from cuppa.runners import PredictionRunner, TrainingRunner, RunnerArgParser
from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH, DEFAULT_CUPPA_CLASSIFIER_PATH
from cuppa.misc.mock_data import MockTrainingData, MockTrainingOutput, MockInputData, MockCvOutput
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary


class TestRunnerArgParser:
    def test_get_args_predict_are_valid(self):

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

    def test_get_args_train_are_valid(self):
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

    def test_run_using_new_input_format(self):

        output_dir = os.path.join(tempfile.gettempdir(), "pycuppa_prediction_run_test")
        os.makedirs(output_dir, exist_ok=True)

        runner = PredictionRunner(
            features_path=MockInputData.path_tsv_new_format_colo,
            sample_id="TEST_SAMPLE",
            classifier_path=DEFAULT_CUPPA_CLASSIFIER_PATH,
            output_dir=output_dir,
            using_old_features_format=False
        )

        runner.run()

        assert os.path.exists(runner.plot_path)
        assert os.path.exists(runner.vis_data_path)
        assert os.path.exists(runner.pred_summ_path)

        shutil.rmtree(output_dir)

    def test_run_on_mock_data_without_outputting_files(self):

        runner = PredictionRunner(
            features_path="/PLACEHOLDER",
            output_dir="/PLACEHOLDER",
            classifier_path=MockTrainingOutput.path_cuppa_classifier
        )

        runner.X = MockTrainingData.X
        runner.get_predictions()
        runner.get_pred_summ()

        assert isinstance(runner.predictions, CuppaPrediction)
        assert isinstance(runner.pred_summ, CuppaPredSummary)


class TestTrainingRunner:

    # from cuppa.test.test_runners import TestTrainingRunner
    # self = TestTrainingRunner

    runner = TrainingRunner(
        features_path="PLACEHOLDER",
        metadata_path="PLACEHOLDER",
        output_dir="/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/output/runner_output/",

        fusion_overrides_path=DEFAULT_FUSION_OVERRIDES_PATH,
        n_jobs=1,
        cache_training=False,
        # log_to_file=True
    )

    runner.X = MockTrainingData.X
    runner.y = MockTrainingData.y
    runner.y_split = MockTrainingData.y_split

    def test_cv_on_mock_data_gives_same_performance(self):
        self.runner.cv_fit()
        self.runner.get_cv_predictions()
        self.runner.get_cv_pred_summ()
        self.runner.get_cv_performance()

        def get_recall(perf):
            return perf.set_index(["class","clf_name"])["recall"].round(4).dropna()

        recall_comparison = pd.DataFrame(dict(
            test = get_recall(self.runner.cv_performance),
            expected = get_recall(MockCvOutput.performance)
        ))

        assert recall_comparison["test"].equals(recall_comparison["expected"])

    def test_train_final_model_successful(self):
        self.runner.train_final_model()

    def _run_individual_steps_on_real_data(self):

        input_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/"

        runner = TrainingRunner(
            features_path=input_dir,
            output_dir="/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/06-pip_env",
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
