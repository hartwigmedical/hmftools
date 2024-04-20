from cuppa.runners.args import RunnerArgParser
from cuppa.runners.prediction_runner import PredictionRunner

if __name__ == "__main__":
    kwargs = RunnerArgParser().get_kwargs_predict()
    runner = PredictionRunner(**kwargs)
    runner.run()
    