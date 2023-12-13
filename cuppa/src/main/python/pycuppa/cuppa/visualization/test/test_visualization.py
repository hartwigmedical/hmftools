import os
from subprocess import CalledProcessError

import pandas as pd
import pytest

from cuppa.misc.mock_data import MockCvOutput
from cuppa.visualization.visualization import CuppaVisDataBuilder, CuppaVisPlotter, CuppaVisData


class TestCuppaVisPlotter:

    def test_r_script_error(self):
        plotter = CuppaVisPlotter(
            CuppaVisData(pd.DataFrame()),
            plot_path="/cuppa_vis.png"
        )

        with pytest.raises(CalledProcessError):
            plotter.plot()

    def _test_plot_one_sample_successful(self):
        predictions = MockCvOutput.predictions_for_vis

        builder = CuppaVisDataBuilder(predictions, sample_id=1)
        vis_data = builder.build()

        plot_path = os.path.expanduser("~/Desktop/cuppa_vis.png")
        vis_data_out_path = os.path.expanduser("~/Desktop/cuppa_vis_data.tsv")

        plotter = CuppaVisPlotter(
            vis_data,
            plot_path=plot_path,
            vis_data_path=vis_data_out_path
        )
        plotter.plot()

        assert os.path.exists(plot_path)
        os.remove(plot_path)
        os.remove(vis_data_out_path)

    def _test_plot_one_sample_from_cuppa_prediction_successful(self):
        predictions = MockCvOutput.predictions_for_vis

        plot_path = os.path.expanduser("~/Desktop/cuppa_vis.png")
        predictions.plot(plot_path=plot_path, sample_id=1)

        assert os.path.exists(plot_path)
        os.remove(plot_path)

