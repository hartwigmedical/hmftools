from subprocess import CalledProcessError

import pytest

from cuppa.misc.executors import RscriptExecutor


class TestRscriptExecutor:
    def test_crash_on_error(self):
        executor = RscriptExecutor(
            args=["dummy_rscript"],
            ignore_error=True
        )
        executor.run()

        assert executor.stderrs[-1].startswith("Fatal error: cannot open file 'dummy_rscript'")

        with pytest.raises(CalledProcessError):
            executor.raise_if_error()

    def test_can_run_r_expression(self):
        executor = RscriptExecutor(
            args=["-e", """ "print('test')" """],
            ignore_error=True
        )
        executor.run()
