from cuppa.runners import PredictionRunner, RunnerArgParser

if __name__ == "__main__":
    kwargs = RunnerArgParser().get_kwargs_predict()
    runner = PredictionRunner(**kwargs)
    runner.run()
    