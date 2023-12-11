import pandas as pd

from cuppa.performance.performance_stats import PerformanceStatsBuilder, PerformanceStats
from cuppa.performance.confusion_matrix import ConfusionMatrix
from cuppa.misc.mock_data import MockCvOutput


class TestPerformanceStatsBuilder:
    def test_build_from_pred_summ(self):
        pred_summ = MockCvOutput.pred_summ
        builder = PerformanceStatsBuilder(pred_summ)

        perf = builder.build()
        assert isinstance(perf, PerformanceStats)


class TestPerformanceStats:

    def test_init_from_tsv(self):
        #df = pd.read_table("/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/02-train_entry_point/cv/report/perf.tsv")
        df = pd.read_table(MockCvOutput.path_performance)
        performance = PerformanceStats.from_data_frame(df)
        assert isinstance(performance, PerformanceStats)

    def test_cuppa_prediction_format_has_correct_indexes(self):
        perf = MockCvOutput.performance
        perf_formatted = perf.to_cuppa_prediction_format()

        expected_index_names = ["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value"]
        assert pd.Series(perf_formatted.index.names).isin(expected_index_names).all()


class TestConfusionMatrix:
    def test_init_from_pred_summ(self):
        pred_summ = MockCvOutput.pred_summ
        confusion = ConfusionMatrix(pred_summ, clf_name="dna_combined")

        assert confusion.counts_matrix.index.name == "pred_class_1" and confusion.counts_matrix.columns.name == "actual_class"
        assert confusion.props_matrix.index.name == "pred_class_1" and confusion.props_matrix.columns.name == "actual_class"

    def _test_plot_succeeds(self):
        pred_summ = MockCvOutput.pred_summ
        confusion = ConfusionMatrix(pred_summ, clf_name="dna_combined")
        confusion.plot()
