from cuppa.runners import RunnerArgParser, TrainingRunner

if __name__ == "__main__":
    kwargs = RunnerArgParser().get_kwargs_train()
    runner = TrainingRunner(**kwargs)
    runner.run()