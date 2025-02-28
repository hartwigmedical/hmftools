import os
import tempfile
import pandas as pd

from tests.mock_data import MockCvOutput, MockVisData
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
            cv_performance=MockCvOutput.performance,
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

    def test_can_plot_one_sample(self):
        plot_path = os.path.join(tempfile.gettempdir(), "cuppa_vis.one_sample.png")

        plotter = CuppaVisPlotter(
            vis_data=MockVisData.vis_data,
            plot_path=plot_path,
        )
        plotter.plot()

        assert os.path.exists(plot_path)
        os.remove(plot_path)

    def test_can_plot_multiple_samples(self):
        plot_path = os.path.join(tempfile.gettempdir(), "cuppa_vis.multi_sample.pdf")

        vis_data_copies = []
        for i in range(0, 3):
            vis_data_copy = MockVisData.vis_data.copy()
            vis_data_copy["sample_id"] = vis_data_copy["sample_id"] + "_" + str(i)
            vis_data_copy = vis_data_copy.drop_cv_performance_data()
            vis_data_copies.append(vis_data_copy)

        vis_data_copies = pd.concat(vis_data_copies)
        vis_data_copies = pd.concat([
            vis_data_copies,
            MockVisData.vis_data.query("data_type=='cv_performance'")
        ])

        plotter = CuppaVisPlotter(
            vis_data=vis_data_copies,
            plot_path=plot_path,
        )

        plotter.plot()

        assert os.path.exists(plot_path)
        os.remove(plot_path)
