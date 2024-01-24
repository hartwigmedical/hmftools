import os
import tempfile
from subprocess import CalledProcessError

import pandas as pd
import pytest

from cuppa.tests.mock_data import MockCvOutput
from cuppa.visualization.visualization import CuppaVisDataBuilder, CuppaVisPlotter, CuppaVisData


class TestCuppaVisDataBuilder:
    def can_build_from_mock_predictions(self):
        pass
        ## TODO: build fails due to missing feature types (i.e. trait). Need to update test feature matrix to not exclude these features
        # builder = CuppaVisDataBuilder(MockCvOutput.predictions, sample_id=1)
        # vis_data = builder.build()
        # assert isinstance(vis_data, CuppaVisData)


class TestCuppaVisPlotter:

    def test_r_script_error(self):

        dummy_vis_data = pd.DataFrame(dict(sample_id="test_sample"), index=[0])
        dummy_vis_data = CuppaVisData(dummy_vis_data)

        plotter = CuppaVisPlotter(
            vis_data=dummy_vis_data,
            plot_path="/cuppa_vis.png"
        )

        with pytest.raises(CalledProcessError):
            plotter.plot_in_r(vis_data_path="/invalid_path", plot_path="/invalid_path")

    def test_can_plot_one_sample(self):
        plot_path = os.path.join(tempfile.gettempdir(), "cuppa_vis.png")

        builder = CuppaVisDataBuilder(MockCvOutput.predictions_for_vis, sample_id=1)
        vis_data = builder.build()

        plotter = CuppaVisPlotter(
            vis_data=vis_data,
            plot_path=plot_path
        )
        plotter.plot()

        assert os.path.exists(plot_path)
        os.remove(plot_path)
