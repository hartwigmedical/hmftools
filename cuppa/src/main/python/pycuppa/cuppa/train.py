from cuppa.runners.args import RunnerArgParser
from cuppa.runners.training_runner import TrainingRunner

if __name__ == "__main__":
    kwargs = RunnerArgParser().get_kwargs_train()
    runner = TrainingRunner(**kwargs)
    runner.run()