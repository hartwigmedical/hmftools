class TestLogging:
    def _test_log_output_manually(self):
        from cuppa.logger import reset_logging_basic_config

        reset_logging_basic_config("/Users/lnguyen/Desktop/train.log")

        from cuppa.runners import TrainingRunner
        from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH
        from cuppa.misc.mock_data import MockTrainingData

        runner = TrainingRunner(
            input_dir="PLACEHOLDER",
            output_dir="/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/output/runner_output/",

            fusion_overrides_path=DEFAULT_FUSION_OVERRIDES_PATH,
            n_jobs=1,
            cache_training=False
        )

        runner.X = MockTrainingData.X
        runner.y = MockTrainingData.y
        runner.y_split = MockTrainingData.y_split

        runner.cv_fit()

        ## Check logs manually @ `output_dir` to check if all messages are outputted there