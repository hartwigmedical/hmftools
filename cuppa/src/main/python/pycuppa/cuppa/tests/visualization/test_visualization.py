import os
import tempfile
from subprocess import CalledProcessError

import pandas as pd
import pytest

from cuppa.tests.mock_data import MockCvOutput, MockVisData
from cuppa.visualization.visualization import CuppaVisDataBuilder, CuppaVisPlotter, CuppaVisData


class TestCuppaVisData:

    EXPECTED_COLUMNS = pd.Series((
        "sample_id", "data_type",
        "clf_group", "clf_name",
        "feat_name", "feat_value",
        "cancer_type", "data_value",
        "rank", "rank_group"
    ))

    def test_can_build_from_mock_predictions(self):

        builder = CuppaVisDataBuilder(
            predictions=MockCvOutput.predictions,
            sample_id=1,
            require_all_feat_types=False ## Mock data is missing 'trait' and 'virus' features
        )
        vis_data = builder.build()

        assert self.EXPECTED_COLUMNS.isin(vis_data.columns).all()
        assert isinstance(vis_data, CuppaVisData)

    def test_can_load_from_tsv(self):
        vis_data = CuppaVisData.from_tsv(MockVisData.path_vis_data)

        assert self.EXPECTED_COLUMNS.isin(vis_data.columns).all()
        assert isinstance(vis_data, CuppaVisData)


class TestCuppaVisPlotter:

    def test_r_script_error(self):

        dummy_vis_data = pd.DataFrame(dict(sample_id="test_sample"), index=[0])
        dummy_vis_data = CuppaVisData(dummy_vis_data)

        plotter = CuppaVisPlotter(
            vis_data=dummy_vis_data,
            plot_path="/dummy_path.pdf"
        )

        ## Override tmp vis data path to be a valid one.
        ## This is usually auto generated from `plot_path`, but here it would not result in a valid path
        plotter._tmp_vis_data_path = os.path.join(tempfile.gettempdir(), "vis_data.tsv")

        with pytest.raises(CalledProcessError):
            plotter.plot()

    def test_can_plot_one_sample(self):
        plot_path = os.path.join(tempfile.gettempdir(), "cuppa_vis.png")

        plotter = CuppaVisPlotter(
            vis_data=MockVisData.vis_data,
            plot_path=plot_path,
        )
        plotter.plot()

        assert os.path.exists(plot_path)
        os.remove(plot_path)
