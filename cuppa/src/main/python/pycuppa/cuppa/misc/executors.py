from __future__ import annotations
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
from cuppa.logger import LoggerMixin


class RscriptExecutor(LoggerMixin):
    def __init__(
        self,
        args: str | list[str],
        ignore_error: bool = False,
        verbose: bool = True
    ):
        self.args = args
        self.ignore_error = ignore_error
        self.verbose = verbose

        self.stderrs: list[str] = []
        self.return_code: int = 0

    @property
    def shell_command(self) -> str:
        if isinstance(self.args, str):
            args = self.args
        elif isinstance(self.args, list):
            args = " ".join(self.args)
        else:
            self.logger.error("`args` must be a str or list[str]")
            raise ValueError

        return "Rscript --vanilla " + args

    def raise_if_error(self) -> None:
        if self.return_code != 0:
            raise CalledProcessError(self.return_code, self.shell_command)

    def run(self) -> None:

        if self.verbose:
            self.logger.info("Running shell command: '%s'" % self.shell_command)

        with Popen(self.shell_command, shell=True, stdout=PIPE, stderr=STDOUT) as process:
            for line in iter(process.stdout.readline, b''):
                stderr = line.decode("utf-8").strip()
                self.logger.info("[R process] " + stderr)

                self.stderrs.append(stderr)

        self.return_code = process.poll()

        if not self.ignore_error:
            self.raise_if_error()