import pandas as pd
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.misc.mock_data import MockProbsPreCal


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

    def test_one_cal_variable_gaussian_kernel(self):
        x = MockProbsPreCal.X["Lung: Small cell"]
        y = MockProbsPreCal.y == "Lung: Small cell"

        transformer = RollingAvgCalibration()
        calibrator = transformer._fit_one_calibrator(
            x, y,
            window_size="variable", n_true_exponent=0.7, edge_weight=0.16, kernel="gaussian"
        )

        assert calibrator.transform([0.3]).round(3) == 0.841

    def test_all_cal_variable_gaussian_kernel(self):
        X = MockProbsPreCal.X
        y = MockProbsPreCal.y

        transformer = RollingAvgCalibration(window_size="variable", n_true_exponent=0.7, edge_weight=0.16, kernel="gaussian")
        transformer.fit(X, y)

        X_new = pd.DataFrame(
            [[0.05, 0.30, 0.65]],
            columns=X.columns
        )

        X_new_trans = transformer.predict_proba(X_new, normalize=False).round(3)
        assert X_new_trans.loc[0, "Lung: Small cell"] == 0.841

        X_new_trans = transformer.predict_proba(X_new, normalize=True).round(3)
        assert X_new_trans.loc[0, "Lung: Small cell"] == 0.297
