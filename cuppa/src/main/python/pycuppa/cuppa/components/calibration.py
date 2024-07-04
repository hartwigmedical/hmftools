from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Iterable, Callable

from numpy._typing import NDArray
from sklearn.base import BaseEstimator
from sklearn.isotonic import IsotonicRegression
from cuppa.components.passthrough import PassthroughRegression
from cuppa.misc.utils import max_col
from cuppa.logger import LoggerMixin


class RollingAvgCalibration(BaseEstimator, LoggerMixin):
    """
    Modified probability calibration approach based on isotonic regression

    Given:
    - an array of probabilities
    - an array of booleans, indicating whether a sample belongs to a cancer type

    Perform the following steps:
    - Sort both arrays by the probabilities
    - Slide a window along the boolean array and take the mean of the values in the window at each step
    - Fit an isotonic regression, with the x-values being the probabilities and y-values being the sliding window means

    This approach results in more granular y-values, which improves the isotonic regression fit

    Parameters
    ----------
    min_true_samples: int
        Minimum number of samples of the target cancer type required to perform the
        calibration. Cancer types with fewer samples than this value will simply have the uncalibrated probabilities
        returned

    window_size: int or 'variable'
        Size of the window (number of samples)

    kernel: 'uniform' or 'gaussian'
        'uniform': mean of the sliding window. 'gaussian': weighted mean where the center has more weight than the
        edges

    edge_weight: float
        Only when kernel='gaussian'. Higher values place more weight on the edges

    n_true_exponent: float
        Only when window_size='variable'. Automatic calculation of the window size based on the formula:
        int(round(n_true_samples ** n_true_exponent))

        Higher values result in wider windows

    save_coords: bool
        Save the probabilities and the rolling average values within each `calibrator` object in
        `self.calibrators`

    bypass: bool
        If True, raw probabilities will be returned upon prediction time (i.e. calibration will be bypassed)

    """

    def __init__(
        self,
        min_true_samples: int = None,
        window_size: int | str = 50,
        kernel: str = "gaussian",
        edge_weight: float = 0.16,
        n_true_exponent: float = 0.7,
        save_coords: bool = False,
        bypass: bool = False
    ):
        self.min_true_samples = min_true_samples
        self.window_size = window_size
        self.edge_weight = edge_weight
        self.n_true_exponent = n_true_exponent
        self.kernel = kernel
        self.save_coords = save_coords
        self.bypass = bypass

    ## Kernel ================================
    def _auto_calc_window_size(self, n_true_samples: int, n_true_exponent: float) -> int:
        return int(round(n_true_samples ** n_true_exponent))

    def _gaussian_kernel(self, window_size=50, edge_weight=0.16) -> NDArray:

        def _normal_dist(x, mean, sd):
            ## Norm dist formula from
            ## https://www.askpython.com/python/normal-distribution
            return (np.pi * sd) * np.exp(-0.5 * ((x - mean) / sd) ** 2)

        weights = _normal_dist(
            x = np.array(range(1, window_size + 1)),
            mean = window_size / 2,
            sd = window_size * edge_weight ## Higher edge_weight values put more weight on indexes near the edge
        )
        weights = weights / np.sum(weights)

        return weights

    def _get_kernel_func(self, kernel, window_size: int = None, edge_weight: float = None) -> Callable:
        if kernel == "uniform":
            func = np.mean

        elif kernel == "gaussian":
            window_weights = self._gaussian_kernel(window_size=window_size, edge_weight=edge_weight)
            func = lambda x: np.sum(x * window_weights)

        else:
            self.logger.error("`kernel` must be 'gaussian' or 'uniform'")
            raise ValueError

        return func

    ## Main ================================
    def _fit_one_calibrator(
        self,
        probs: Iterable[float],
        bools: Iterable[bool],
        min_true_samples: int = None,
        window_size: str | int = "variable",
        n_true_exponent: float = 0.7,
        edge_weight: float = 0.16,
        kernel: str = "gaussian",
        save_coords: bool = False
    ) -> IsotonicRegression | PassthroughRegression:

        ## Checks --------------------------------
        probs = np.array(probs)
        bools = np.array(bools)
        n_true_samples = sum(bools)

        if probs.ndim != 1:
            self.logger.error("`probs` must be an 1d array")
            raise TypeError

        if bools.ndim != 1 or bools.dtype != "bool":
            self.logger.error("`bool` must be an 1d boolean array")
            raise TypeError

        if edge_weight <= 0:
            self.logger.error("`edge_weight` must be >0")
            raise ValueError

        if not isinstance(window_size, str) and window_size <= 0:
            self.logger.error("`window_size` must be >0")
            raise ValueError

        ## Prepare kernel --------------------------------
        if window_size == "variable":
            window_size = self._auto_calc_window_size(n_true_samples=n_true_samples, n_true_exponent=n_true_exponent)

        kernel_func = self._get_kernel_func(kernel=kernel, window_size=window_size, edge_weight=edge_weight)

        ## Main --------------------------------
        ## Sort data by prob
        coords = pd.DataFrame(dict(bool=bools, prob=probs), index=None)
        coords.sort_values("prob", ascending=False, inplace=True)

        ## Ignore calibration if there are not enough samples
        if min_true_samples is not None and n_true_samples < min_true_samples:
            calibrator = PassthroughRegression()

        else:
            ## Add padding to bool array to handle kernel going beyond edges
            bool_padded = np.concatenate([
                np.ones(window_size),
                coords["bool"],
                np.zeros(window_size)
            ])

            ## Rolling averages to smooth out bool array
            avgs = pd.Series(bool_padded).rolling(window=window_size, center=True).apply(kernel_func)

            ## Discard padding
            n_samples = len(bools)
            avgs = avgs[window_size:(window_size + n_samples)]

            coords["avg"] = avgs.values

            calibrator = IsotonicRegression(out_of_bounds="clip")
            calibrator.fit(coords["prob"], coords["avg"])

        ## Store training x,y coords for debugging
        if save_coords:
            calibrator.coords = coords

        return calibrator

    def fit(self, X: pd.DataFrame, y: pd.Series) -> "RollingAvgCalibration":
        """

        Parameters
        ----------
        X: pandas DataFrame
            Probabilities for each sample (row) and cancer type (columns)

        y: pandas Series of type str or Categorical
            Cancer type labels for each sample

        Returns
        -------
        self
            fitted estimator
        """
        if not isinstance(X, pd.DataFrame):
            self.logger.error("`X` must be a pandas dataframe")
            raise ValueError

        calibrators = {}
        for class_i in np.unique(y):
            calibrator_i = self._fit_one_calibrator(
                X[class_i],
                y==class_i,
                min_true_samples = self.min_true_samples,
                window_size = self.window_size,
                edge_weight = self.edge_weight,
                n_true_exponent = self.n_true_exponent,
                kernel = self.kernel,
                save_coords = self.save_coords
            )
            calibrators[class_i] = calibrator_i

        self.calibrators = calibrators

        return self

    def predict_proba(self, X: pd.DataFrame, normalize: bool = True) -> pd.DataFrame:

        if self.bypass or X.shape[0] == 0:
            return X

        if not isinstance(X, pd.DataFrame):
            self.logger.error("`X` must be a pandas dataframe")
            raise ValueError

        probs_cal = {}
        for class_name in X.columns:
            calibrator = self.calibrators[class_name]
            probs_cal[class_name] = calibrator.predict(X[class_name])

        probs_cal = pd.DataFrame(probs_cal, index=X.index)

        if normalize:
            ## Scale probs from 0 to 1
            probs_cal = probs_cal.T / np.sum(probs_cal, axis=1)
            probs_cal = probs_cal.T

        return probs_cal

    def predict(self, X: pd.DataFrame, y: pd.Series = None) -> pd.Series:
        return max_col(self.predict_proba(X))

    def transform(self, X: pd.DataFrame, y: pd.Series = None) -> pd.DataFrame:
        return self.predict_proba(X)

    def fit_transform(self, X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
        return self.fit(X, y).transform(X)

    def set_output(self, transform: str = None) -> "RollingAvgCalibration":
        return self

    # Calibration curves ================================
    def get_cal_curves(self, x_resolution=500, long_format=True) -> pd.DataFrame:

        ## Get y values
        x_values = np.linspace(0, 1, x_resolution + 1)

        cal_curves = {
            class_name: calibrator.predict(x_values)
            for class_name, calibrator in self.calibrators.items()
        }

        ## Wide format with x values as first col
        cal_curves = pd.DataFrame(cal_curves)
        cal_curves["x"] = x_values

        if long_format:
            cal_curves = cal_curves.melt(id_vars="x", var_name="class", value_name="y")
            cal_curves = cal_curves[['class', 'x', 'y']]

        return cal_curves

    ## TODO: Make plot using ggplot2 in R