import sys

from cuppa.runners.args import RunnerArgParser
from cuppa.runners.prediction_runner import PredictionRunner
from cuppa.runners.training_runner import TrainingRunner


class TestRunnerArgParser:
    def test_can_add_and_parse_one_arg(self):

        sys.argv = [
            "example_module.py",
            "--features_path", "/path/to/features/",
        ]

        arg_parser = RunnerArgParser()
        arg_parser.add_argument("features_path")

        kwargs = arg_parser.get_kwargs()

        assert kwargs["features_path"] == "/path/to/features/"

    def test_can_create_prediction_runner_from_required_command_line_args(self):
        sys.argv = [
            "example_module.py",
            "--features_path", "/path/to/features/",
            "--sample_id", "SAMPLE",
            "--output_dir", "/path/to/output/",
            "--classifier_path", "/path/to/cuppa_classifier.pickle.gz",
            "--cv_predictions_path", "/path/to/cv_predictions.tsv.gz",
        ]

        arg_parser = RunnerArgParser()
        kwargs = arg_parser.get_kwargs_predict()
        runner = PredictionRunner(**kwargs)
        assert True

    def test_can_create_training_runner_from_required_command_line_args(self):
        sys.argv = [
            "example_module.py",
            "--features_path", "/path/to/features/",
            "--metadata_path", "/path/to/metadata_file",
            "--output_dir", "/path/to/output/",
        ]

        arg_parser = RunnerArgParser()
        kwargs = arg_parser.get_kwargs_train()

        runner = TrainingRunner(**kwargs)
        assert True
