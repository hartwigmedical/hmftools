import pandas as pd
import os
import tempfile

from tests.mock_data import MockCvOutput
from cuppa.performance.performance_stats import PerformanceStatsBuilder, PerformanceStats
from cuppa.performance.confusion_matrix import ConfusionMatrix


class TestPerformanceStatsBuilder:
    def test_can_build_from_pred_summ(self):
        pred_summ = MockCvOutput.pred_summ
        builder = PerformanceStatsBuilder(pred_summ)

        perf = builder.build()

        assert list(perf.columns) == ['class', 'clf_name', 'n_total', 'n_predicted', 'n_correct', 'recall', 'precision']
        assert isinstance(perf, PerformanceStats)


class TestPerformanceStats:

    def test_can_load_from_tsv(self):
        performance = PerformanceStats.from_tsv(MockCvOutput.path_performance)
        assert isinstance(performance, PerformanceStats)


class TestConfusionMatrix:
    def test_can_initialize_from_prediction_summary(self):
        pred_summ = MockCvOutput.pred_summ
        confusion = ConfusionMatrix(pred_summ, clf_name="dna_combined")

        assert confusion.counts_matrix.index.name == "pred_class_1" and confusion.counts_matrix.columns.name == "actual_class"
        assert confusion.props_matrix.index.name == "pred_class_1" and confusion.props_matrix.columns.name == "actual_class"

    def test_can_plot_confusion_matrix(self):
        pred_summ = MockCvOutput.pred_summ
        confusion = ConfusionMatrix(pred_summ, clf_name="dna_combined")

        plot_path = os.path.join(tempfile.gettempdir(), "plot.pdf")
        confusion.plot(plot_path)

        assert os.path.exists(plot_path)

        os.remove(plot_path)
