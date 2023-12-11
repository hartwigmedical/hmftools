from functools import cached_property
from typing import Literal

import numpy as np
import pandas as pd
from numpy._typing import NDArray
from scipy.stats import false_discovery_control
from sklearn.base import BaseEstimator
from sklearn.feature_selection import chi2

from cuppa.components.preprocessing import MaxScaler
from cuppa.logger import LoggerMixin


class Chi2FeatureSelector(BaseEstimator, LoggerMixin):

    def __init__(
        self,
        mode: Literal["pvalue", "qvalue", "fdr", "rank", "k_best"] = "fdr",
        threshold: float | int = 0.001,
    ):
        """
        Chi2 feature selection using the chi2() sklearn function

        The chi2 test quantifies how much each feature varies across every class (i.e. cancer type) in the training set
        compared to an expected (background) value for that feature. It is an extremely fast multi-class feature
        selection method compared to one-vs-rest pairwise testing (e.g. using wilcoxon rank sum test).

        Parameters
        ----------
        mode: str
            The type of measure to use when applying the `threshold` value. Can be:
            pvalue: p-value
            qvalue: p-value adjusted for multiple testing using the Benjamini-Hochberg method
            fdr: Same as q-value
            rank: Top features
            k_best: Same as rank

        threshold:
            The cutoff to use for the measure specified in `mode`

        """
        self.mode = mode
        self.threshold = threshold

        self._check_mode_valid()
        self._check_threshold_valid()

        self.test_results = None

    def _check_mode_valid(self):
        valid_modes = [
            "pvalue",
            "qvalue","fdr",
            "rank","k_best"
        ]

        if self.mode not in valid_modes:
            self.logger.error("Valid modes are: " + ", ".join(valid_modes))
            raise ValueError

    def _check_threshold_valid(self):
        if self.threshold < 0:
            self.logger.error("`threshold` must be non-negative")
            raise ValueError

        if self.mode in ["rank", "k_best"] and not isinstance(self.threshold, int):
            self.logger.error("`threshold` must be an int when `mode` is 'rank' or 'k_best'")
            raise ValueError

    def fit(self, X: pd.DataFrame, y: pd.Series):

        stat_values, pvalues = chi2(X, y)
        pvalues[np.isnan(pvalues)] = 1
        qvalues = false_discovery_control(pvalues, method="bh")

        results = pd.DataFrame(
            dict(
                stat = stat_values,
                pvalue = pvalues,
                qvalue = qvalues
            ),
            index = X.columns
        )

        results["rank"] = results["stat"].rank(method="first", ascending=False, na_option="bottom").astype(int)

        self.test_results = results

    @cached_property
    def selected_features(self) -> NDArray:

        if self.test_results is None:
            self.logger.error(self.__class__.__name__ + "is not yet fitted")
            raise LookupError

        feature_names = self.test_results.index

        if self.mode == "pvalue":
            target_values = self.test_results["pvalue"]
            selected_features = feature_names[target_values < self.threshold]

        if self.mode in ["qvalue", "fdr"]:
            target_values = self.test_results["qvalue"]
            selected_features = feature_names[target_values < self.threshold]

        if self.mode in ["rank", "k_best"]:
            target_values = self.test_results["rank"]
            selected_features = feature_names[target_values <= self.threshold]

        return selected_features

    def transform(self, X: pd.DataFrame, y = None):
        return X[self.selected_features]

    def fit_transform(self, X: pd.DataFrame, y = None):
        self.fit(X, y)
        return self.transform(X)

    def set_output(self, transform = None):
        return self


class MaxScaledChi2FeatureSelector(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        mode: Literal["pvalue", "qvalue", "fdr", "rank", "k_best"] = "fdr",
        threshold: float | int = 0.001,
        clip: bool = True
    ):
        """
        MaxScaler followed by Chi2FeatureSelector

        This class exists as the division operation in MaxScaler is computationally expensive when dealing with 50,000+
        features, as is the case for the gene expression and alt SJ features.

        When calling `transform()`, the input dataframe is subsetted for the select features. Then max scaling is
        performed only on these features. This avoids max scaling from being unnecessarily performed on unused features.

        Parameters
        ----------
        mode: str
            The type of measure to use when applying the `threshold` value. Can be:
            pvalue: p-value
            qvalue: p-value adjusted for multiple testing using the Benjamini-Hochberg method
            fdr: Same as q-value
            rank: Top features
            k_best: Same as rank

        threshold:
            The cutoff to use for the measure specified in `mode`

        clip : bool
            Set values <0 and >1 after division to 0 and 1 respectively during max scaling?
        """

        self.threshold = threshold
        self.mode = mode
        self.clip = clip

        self.max_scaler = MaxScaler(clip=clip)
        self.chi2_selector = Chi2FeatureSelector(threshold=threshold, mode=mode)

        self.selected_features = None

    def fit(self, X: pd.DataFrame, y: pd.Series):
        X_scaled = self.max_scaler.fit_transform(X, y)
        X_selected = self.chi2_selector.fit_transform(X_scaled, y)

        self.selected_features = self.chi2_selector.selected_features
        self.max_scaler.max_values = self.max_scaler.max_values[self.selected_features]

    def transform(self, X: pd.DataFrame, y: pd.Series):
        X_selected = X[self.selected_features]
        X_scaled = self.max_scaler.transform(X_selected)
        return X_scaled

    def fit_transform(self, X: pd.DataFrame, y = None):
        self.fit(X, y)
        return self.transform(X)

    def set_output(self, transform = None):
        return self
