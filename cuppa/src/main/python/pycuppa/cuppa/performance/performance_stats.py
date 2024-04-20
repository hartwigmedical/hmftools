from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy._typing import NDArray

from cuppa.constants import PAN_CANCER_CLASS_NAME
from cuppa.logger import LoggerMixin
from cuppa.misc.utils import check_required_columns, bin_array

if TYPE_CHECKING:
    from cuppa.classifier.cuppa_prediction import CuppaPredSummary


class PerformanceStatsBuilder(LoggerMixin):
    def __init__(
        self,
        pred_summ: CuppaPredSummary,
        prob_bins: NDArray = np.linspace(0, 1, 6)
    ):
        """
        Creates a PerformanceStats object from a CuppaPredSummary object

        Parameters
        ----------
        pred_summ: CuppaPredSummary, child of pandas dataframe
            Dataframe summarizing the top-N predictions

        prob_bins: array-like
            Specifies the boundaries for splitting probabilities by bin, e.g.: [0.0 , 0.2, 0.4, 0.6, 0.8, 1.0].
            Only used when calling `build(by_prob_bin=True)`

        """

        self.pred_summ = pred_summ
        self.prob_bins = prob_bins

    def _check_actual_class_column_exists(self) -> None:
        if not self.pred_summ._has_actual_class_column:
            self.logger.error("Cannot calculate performance when `actual_class` column is absent")
            raise KeyError

    @property
    def clf_names(self) -> NDArray:
        return self.pred_summ.clf_names

    @property
    def classes(self) -> NDArray:
        return self.pred_summ.classes

    @cached_property
    def _pred_summ_parsed(self) -> pd.DataFrame:

        self._check_actual_class_column_exists()

        ## Select relevant columns
        df = self.pred_summ[["clf_name", "actual_class", "pred_class_1", "pred_prob_1", "is_correct_pred"]]

        ## Remove NA predictions
        df = df[~df["pred_class_1"].isna()]

        ## Set classes
        df["actual_class"] = pd.Categorical(df["actual_class"], self.classes)
        df["pred_class_1"] = pd.Categorical(df["pred_class_1"], self.classes)

        return df

    @staticmethod
    def _get_stats_per_class(_pred_summ_parsed: pd.DataFrame) -> pd.DataFrame:

        stats = {}

        stats["n_total"] = _pred_summ_parsed.groupby("clf_name", observed=False)["actual_class"].value_counts()
        stats["n_predicted"] = _pred_summ_parsed.groupby("clf_name", observed=False)["pred_class_1"].value_counts()
        stats["n_correct"] = _pred_summ_parsed.groupby(["clf_name","actual_class"], observed=False)["is_correct_pred"].sum()

        ## Combine into one dataframe
        for stats_i in stats.values():
            stats_i.index.names = ["clf_name", "class"]
        stats = pd.DataFrame(stats)

        stats = stats.reset_index()

        return stats

    @staticmethod
    def _get_stats_pan_class(_pred_summ_parsed: pd.DataFrame) -> pd.DataFrame:

        stats = _pred_summ_parsed \
            .groupby(["clf_name"]) \
            .agg(dict(actual_class="count", is_correct_pred="sum"))

        stats.columns = ["n_total", "n_correct"]

        stats["n_predicted"] = np.nan ## Number predicted pan-cancer is undefined
        stats["class"] = PAN_CANCER_CLASS_NAME

        stats.index.name = "clf_name"
        stats = stats.reset_index()

        return stats

    @staticmethod
    def _add_metrics(stats) -> None:
        ## Convert to float to allow divide by zero
        stats["recall"] = stats["n_correct"].astype(float) / stats["n_total"].astype(float)
        stats["precision"] = stats["n_correct"].astype(float) / stats["n_predicted"].astype(float)

    COLUMN_NAMES = ["class", "clf_name", "n_total", "n_predicted", "n_correct", "recall", "precision"]

    def _get_performance(self, _pred_summ_parsed: pd.DataFrame) -> pd.DataFrame:

        stats = pd.concat([
            self._get_stats_pan_class(_pred_summ_parsed),
            self._get_stats_per_class(_pred_summ_parsed)
        ])

        self._add_metrics(stats)
        stats["n_predicted"] = stats["n_predicted"].fillna(0).astype(int)

        ## Rearrange columns
        stats = stats[self.COLUMN_NAMES]

        ## Force classifier order
        stats = stats.sort_values(["class", "clf_name"])

        return stats

    def _get_performance_by_prob_bin(self, _pred_summ_parsed) -> PerformanceStats:

        df = _pred_summ_parsed

        prob_bins = bin_array(df["pred_prob_1"], self.prob_bins)
        df["prob_bin"] = prob_bins

        stats = []
        for prob_bin in prob_bins.categories:
            #prob_bin = "(0.6,0.8]"
            df_i = df[df["prob_bin"] == prob_bin]
            stats_i = self._get_performance(df_i)
            stats_i.insert(stats_i.columns.get_loc("clf_name") + 1, "prob_bin", prob_bin)
            stats.append(stats_i)

        stats = pd.concat(stats)

        stats["prob_bin"] = pd.Categorical(stats["prob_bin"], prob_bins.categories)

        stats = stats.sort_values("prob_bin").sort_values(["class", "clf_name"])
        stats = stats.reset_index(drop=True)

        return PerformanceStats.from_data_frame(stats)

    def build(self, by_prob_bin: bool = False) -> PerformanceStats:
        if not by_prob_bin:
            df = self._get_performance(self._pred_summ_parsed)
        else:
            df = self._get_performance_by_prob_bin(self._pred_summ_parsed)

        return PerformanceStats(df)


class PerformanceStats(pd.DataFrame, LoggerMixin):
    def __init__(self, performance: pd.DataFrame, *args, **kwargs):
        """
        A table showing performance stats

        It has the following columns:
        - class: predicted/actual class
        - clf_name: classifier name
        - n_total: total number of samples per class
        - n_predicted: number of samples predicted as the respective class
        - n_correct: number of samples correctly predicted as the respective class
        - recall: fraction of samples correctly predicted as the respective class
        - precision: of the samples predicted as a certain class, the fraction correctly predicted

        This class extends pandas dataframe with additional methods.
        - Please use the methods from_data_frame() or from_tsv() to instantiate this class.
        - For exporting, use the to_tsv() method.
        - The to_cuppa_prediction_format() method is only used internally for building a CuppaPrediction object.

        Parameters
        ----------
        performance: pandas dataframe
            Table with the columns described above
        """
        super().__init__(performance, *args, **kwargs)

    @property
    def _constructor(self):
        return PerformanceStats

    @property
    def is_by_prob_bin(self) -> bool:
        return "prob_bin" in self.columns

    @property
    def classes(self) -> pd.Series:
        return self["class"].unique()

    @property
    def clf_names(self) -> pd.Series:
        return self["clf_name"].unique()

    def to_tsv(self, path: str, verbose: bool = False) -> None:
        if verbose:
            self.logger.info("Writing performance to: " + path)

        self.to_csv(path, sep="\t", index=False)

    @staticmethod
    def from_data_frame(df: pd.DataFrame) -> PerformanceStats:
        required_columns = ["class", "clf_name", "n_total", "n_predicted", "n_correct", "recall", "precision"]
        check_required_columns(df, required_columns=required_columns)
        return PerformanceStats(df)

    @classmethod
    def from_tsv(cls, path: str) -> PerformanceStats:
        df = pd.read_table(path)
        return cls.from_data_frame(df)
