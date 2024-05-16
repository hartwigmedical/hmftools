from __future__ import annotations

import os
import tempfile
from functools import cached_property
from typing import Optional

import numpy as np
import pandas as pd

from cuppa.classifier.cuppa_prediction import CuppaPredSummary
from cuppa.constants import PAN_CANCER_CLASS_NAME, RSCRIPT_PLOT_CONFUSION_PATH
from cuppa.logger import LoggerMixin
from cuppa.misc.executors import RscriptExecutor
from cuppa.performance.performance_stats import PerformanceStats


class ConfusionMatrix(LoggerMixin):
    def __init__(self, pred_summ: CuppaPredSummary, clf_name: str = None):
        """
        Confusion matrix

        A matrix where rows are predicted classes and columns are actual classes. The diagonal shows per actual class
        the correct classifications, and the off-diagonal shows the incorrect classifications.

        Parameters
        ----------
        pred_summ: CuppaPredSummary, child of pandas DataFrame
            Dataframe summarizing the top-N predictions

        clf_name: str
            Classifier name. When `pred_summ` contains the predictions from multiple classifiers, `clf_name` needs to be
            specified as a confusion matrix can only made per classifier

        Properties
        ----------
            counts_matrix: Number of correct/incorrect classifications per actual class
            props_matrix: Proportion of correct/incorrect classifications per actual class


        Methods
        -------
            plot: Plot confusion matrix as a heatmap correct/incorrect proportion as the fill colors, and the
            correct/incorrect counts as the labels.

        """

        self.clf_name = clf_name
        self.pred_summ = pred_summ

        self._select_one_clf()

    def _select_one_clf(self) -> None:
        n_clfs = self.pred_summ.clf_names.__len__()

        if n_clfs == 1:
            self.clf_name = self.pred_summ["clf_name"].iloc[0]
            return None

        if self.clf_name is None and n_clfs > 1:
            self.logger.error("`pred_summ` contains results of multiple classifiers. Please specify `clf_name`")
            raise Exception

        if self.clf_name not in self.pred_summ.clf_names:
            self.logger.error("Valid values for `clf_name` are: " + ", ".join(self.pred_summ.clf_names))
            raise ValueError

        self.pred_summ = self.pred_summ[self.pred_summ["clf_name"] == self.clf_name]


    @cached_property
    def _performance(self) -> PerformanceStats:
        perf = self.pred_summ.performance()
        perf = perf[perf["clf_name"] == self.clf_name]
        return perf

    def _get_stat_all(self, column: str) -> pd.DataFrame:
        stat = self._performance.loc[self._performance["class"] == PAN_CANCER_CLASS_NAME, column]
        return stat.squeeze()

    def _make_confusion_matrix(self, normalize: bool = False) -> pd.DataFrame:
        ## Subset input data
        df = self.pred_summ[["actual_class","pred_class_1"]]

        ## Remove rows where sample does not have a prediction
        df = df[~df["pred_class_1"].isna()]

        ## Get unique classes
        classes = []
        classes += df["actual_class"].unique().tolist()
        classes += df["pred_class_1"].unique().tolist()
        classes = pd.Series(classes).sort_values().unique()

        df["actual_class"] = pd.Categorical(df["actual_class"], classes)
        df["pred_class_1"] = pd.Categorical(df["pred_class_1"], classes)

        ## Confusion matrix
        matrix = pd.crosstab(
            df["pred_class_1"], df["actual_class"],
            normalize=False if not normalize else "columns"
        )

        ## Add pan-cancer stat
        matrix.insert(0, PAN_CANCER_CLASS_NAME, np.nan) ## Add empty rows

        matrix = matrix.T ## Add empty columns
        matrix.insert(0, PAN_CANCER_CLASS_NAME, np.nan)
        matrix = matrix.T

        stat_all = self._get_stat_all("n_correct")
        if normalize:
            stat_all /= self._get_stat_all("n_total")
        matrix.iloc[0, 0] = stat_all

        if not normalize:
            matrix = matrix.astype("Int64")

        return matrix

    @cached_property
    def counts_matrix(self) -> pd.DataFrame:
        return self._make_confusion_matrix(normalize=False)

    @cached_property
    def props_matrix(self) -> pd.DataFrame:
        return self._make_confusion_matrix(normalize=True)

    def __repr__(self) -> str:
        props_and_counts = self.counts_matrix.astype(str) + " (" + self.props_matrix.round(2).astype(str) + ")"
        return props_and_counts.__repr__()

    def plot(
        self,
        plot_path: str,
        width: int = None, ## inches
        height: int = None, ## inches
    ) -> None:

        tmp_pred_summ_path = os.path.join(os.path.dirname(plot_path), "cuppa.pred_summ.tmp.tsv.gz")

        self.logger.debug("Writing pred summ to temporary path: " + tmp_pred_summ_path)
        self.pred_summ.to_tsv(tmp_pred_summ_path)

        try:
            executor = RscriptExecutor([
                RSCRIPT_PLOT_CONFUSION_PATH,
                tmp_pred_summ_path,
                self.clf_name,
                plot_path,
                str(width) if width is not None else "",
                str(height) if height is not None else ""
            ])
            executor.run()
        finally:
            os.remove(tmp_pred_summ_path)
