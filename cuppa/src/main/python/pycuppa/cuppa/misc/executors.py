from __future__ import annotations

import subprocess
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
from cuppa.logger import LoggerMixin


class RscriptExecutor(LoggerMixin):
    def __init__(
        self,
        args: list[str],
        ignore_error: bool = False, ## Used for testing
        verbose: bool = True
    ):
        self.args = args
        self.ignore_error = ignore_error
        self.verbose = verbose

        self.stderrs: list[str] = []
        self.return_code: int = 0

        self.check_R_installed()

    def check_R_installed(self) -> None:

        result = subprocess.run(["which", "Rscript"], capture_output=True)

        if result.returncode > 0:
            self.logger.error("R is not installed!")
            raise CalledProcessError
        else:
            self.logger.debug("Found Rscript binary at: " + result.stdout.decode("utf-8").strip())

    @property
    def command(self) -> list[str]:

        if not isinstance(self.args, list):
            self.logger.error("`args` must be a list[str]")
            raise ValueError

        return ["Rscript", "--vanilla"] + self.args

    def run(self) -> None:

        if self.verbose:
            self.logger.info("Running command: " + " ".join(self.command))

        with Popen(args=self.command, shell=False, stdout=PIPE, stderr=STDOUT) as process:
            for line in iter(process.stdout.readline, b''):
                stderr = line.decode("utf-8").strip()
                self.logger.debug("[R] " + stderr)

                self.stderrs.append(stderr)

        self.return_code = process.poll()

        if not self.ignore_error:
            self.raise_if_error()

    def raise_if_error(self) -> None:
        if self.return_code != 0:
            raise CalledProcessError(self.return_code, " ".join(self.command))
