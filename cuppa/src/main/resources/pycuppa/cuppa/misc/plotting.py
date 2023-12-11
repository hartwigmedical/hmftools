import os
import matplotlib
import matplotlib.pyplot as plt
from typing import Optional

from cuppa.logger import LoggerMixin

class PlotnineFigExporter(LoggerMixin):
    def __init__(
        self,
        width: int | float = 16,
        height: int | float = 10,
        dpi: int = 300,
        verbose: bool = True
    ):
        self.width = width
        self.height = height
        self.dpi = dpi
        self.verbose = verbose

    def _check_input(self, fig: matplotlib.figure.Figure):

        if isinstance(fig, matplotlib.figure.Figure):
            return fig

        self.logger.error("`fig` must be of class `matplotlib.figure.Figure`")
        raise TypeError

    def _check_path(self, path: str):
        if path is None:
            path = os.path.expanduser('~/Desktop/plot.pdf')
            self.logger.warning("No path was specified. Using: " + path)

        return path

    def export(
        self,
        fig: matplotlib.figure.Figure,
        path: Optional[str]
    ):
        path = self._check_path(path)

        fig = self._check_input(fig)
        fig.set_size_inches(self.width, self.height)
        fig.savefig(path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
