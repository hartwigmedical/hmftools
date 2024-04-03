import pandas as pd

from tests.mock_data import MockCvOutput, MockTrainingData
from cuppa.components.calibration import RollingAvgCalibration


class TestKernel:

    def test_edge_weight_default(self):
        kernel = RollingAvgCalibration()._gaussian_kernel(window_size=20, edge_weight=0.16)
        assert round(kernel[0], 4) == 0.0024

    def test_edge_weight_high(self):
        kernel = RollingAvgCalibration()._gaussian_kernel(window_size=20, edge_weight=1)
        assert round(kernel[0], 4) == 0.0471

    def test_auto_calc_window_size_small(self):
        window_size = RollingAvgCalibration()._auto_calc_window_size(n_true_samples=5, n_true_exponent=0.7)
        assert window_size == 3

    def test_auto_calc_window_size_large(self):
        window_size = RollingAvgCalibration()._auto_calc_window_size(n_true_samples=5, n_true_exponent=1.5)
        assert window_size == 11


class TestRollingAvgCalibration:

    def test_one_class_calibration_gives_expected_results(self):
        probs = MockCvOutput.probs_per_clf["gen_pos"]["Lung"]
        y = MockTrainingData.y=="Lung"

        transformer = RollingAvgCalibration()
        calibrator = transformer._fit_one_calibrator(
            probs, y,
            window_size="variable", n_true_exponent=0.7, edge_weight=0.16, kernel="gaussian"
        )

        probs_new = [0.1, 0.2, 0.3, 0.4]
        probs_new_cal = calibrator.transform(probs_new).round(3).tolist()
        assert probs_new_cal == [0.0, 0.348, 0.765, 0.918]

    def test_multi_class_calibration_gives_expected_results(self):
        probs = MockCvOutput.probs_per_clf["gen_pos"]
        y = MockTrainingData.y

        transformer = RollingAvgCalibration(window_size="variable", n_true_exponent=0.7, edge_weight=0.16, kernel="gaussian")
        transformer.fit(probs, y)

        probs_new = pd.DataFrame(
            [[0.05, 0.30, 0.65, 0, 0]],
            columns=probs.columns
        )

        probs_new_cal = transformer\
            .predict_proba(probs_new, normalize=False).round(3) \
            .iloc[0].values.tolist()

        assert probs_new_cal == [0.0, 0.345, 1.0, 0.0, 0.0]
