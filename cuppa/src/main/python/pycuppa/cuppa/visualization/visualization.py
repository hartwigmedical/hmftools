from __future__ import annotations

import os

import numpy as np
import pandas as pd
from functools import cached_property

from cuppa.constants import RSCRIPT_PLOT_PREDICTIONS_PATH, META_CLF_NAMES, CLF_GROUPS
from cuppa.logger import LoggerMixin
from cuppa.misc.executors import RscriptExecutor
from cuppa.misc.utils import check_required_columns

from typing import TYPE_CHECKING, Optional


if TYPE_CHECKING:
    from cuppa.classifier.cuppa_prediction import CuppaPrediction
    from cuppa.performance.performance_stats import PerformanceStats


class CuppaVisDataBuilder(LoggerMixin):
    def __init__(
        self,
        predictions: CuppaPrediction,
        cv_performance: PerformanceStats,
        sample_id: Optional[str | int] = None,
        require_all_feat_types: bool = True,

        min_driver_likelihood: Optional[float] = 0.2, ## LINX driver likelihood
        min_driver_contrib: Optional[float] = 0.2, ## log(odds)

        verbose: bool = True
    ):
        self.predictions = predictions
        self.cv_performance = cv_performance
        self.sample_id = sample_id
        self.require_all_feat_types = require_all_feat_types

        self.min_driver_likelihood = min_driver_likelihood
        self.min_driver_contrib = min_driver_contrib

        self.verbose = verbose

        self._check_snv_counts_feature_exists()

    COLUMN_NAMES = [
        "sample_id", "data_type", "clf_group", "clf_name",
        "feat_name", "feat_value", "cancer_type", "data_value",
        "rank", "rank_group"
    ]

    FEAT_NAME_SNV_COUNT = "event.tmb.snv_count"

    def _check_snv_counts_feature_exists(self):
        if self.FEAT_NAME_SNV_COUNT not in self.predictions.index.get_level_values("feat_name"):
            self.logger.error(
                "Feature contributions / values for '%s' do not exist. "
                "This is required for calculating relative mutational signature counts"
            )
            raise LookupError

    ## Probs ================================
    def get_probs(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting probabilities")

        return self.predictions.get_data_types("prob").wide_to_long()

    ## Signatures ================================
    def get_signatures(self) -> pd.DataFrame:
        if self.verbose:
            self.logger.debug("Getting mutational signatures")

        return self.predictions.get_data_types("sig_quantile").wide_to_long()

    ## Top event features ================================
    @cached_property
    def feat_contrib(self) -> CuppaPrediction:
        return self.predictions.get_data_types("feat_contrib")

    def _subset_feat_contrib_by_feat_pattern(self, pattern: str) -> CuppaPrediction | None:
        feat_contrib = self.feat_contrib[
            self.feat_contrib.index.get_level_values("feat_name").str.contains(pattern)
        ]

        if len(feat_contrib) == 0:
            if self.require_all_feat_types:
                self.logger.error("No features found with regex '%s'" % pattern)
                raise LookupError
            else:
                self.logger.warning("No features found with regex '%s'. Returning None" % pattern)
                return None

        return feat_contrib

    def _get_top_drivers(self) -> CuppaPrediction | None:
        feat_contrib = self._subset_feat_contrib_by_feat_pattern("driver")

        if feat_contrib is None:
            return None

        index = feat_contrib.index.to_frame(index=False)
        row_maxes = (feat_contrib.max(axis=1) >= self.min_driver_contrib).values
        selected_rows = (row_maxes >= self.min_driver_contrib) | (index["feat_value"] >= self.min_driver_likelihood)

        return feat_contrib[selected_rows.values]

    def _get_existing_features(self, feat_pattern: str) -> CuppaPrediction | None:
        feat_contrib = self._subset_feat_contrib_by_feat_pattern(feat_pattern)

        if feat_contrib is None:
            return None

        max_contrib_by_row = feat_contrib.max(axis=1)
        feat_values = feat_contrib.index.get_level_values("feat_value")
        return feat_contrib[(max_contrib_by_row > 0) & (feat_values > 0)]

    def _get_existing_fusions(self) -> CuppaPrediction | None:
        return self._get_existing_features("fusion")

    def _get_existing_viruses(self) -> CuppaPrediction | None:
        return self._get_existing_features("virus")

    def get_top_feat_contrib(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting top event features")

        wide = self.predictions.__class__.concat([  ## Indirect call to CuppaPrediction.concat() to avoid circular import
            self._subset_feat_contrib_by_feat_pattern("tmb"),
            self._subset_feat_contrib_by_feat_pattern("trait"),
            self._get_existing_fusions(),
            self._get_existing_viruses(),
            self._get_top_drivers(),
            self._subset_feat_contrib_by_feat_pattern("sv")
        ])

        long = wide.wide_to_long()

        ## Calculate odds
        long["data_value"] = np.exp(1) ** long["data_value"]

        ## Strip clf_name from feat_name
        long["feat_name"] = long["feat_name"].str.split(".", n=1).map(lambda x: x[1])

        return long

    ## CV performance ================================
    def parse_cv_performance(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Parsing cross-validation performance data")

        perf = self.cv_performance.copy()

        ## Cast to long format
        perf = perf.melt(id_vars=["class", "clf_name"], var_name="feat_name")

        perf = perf.rename(columns={
            "class": "cancer_type",
            "value": "data_value"
        })

        ## Subset for relevant rows
        perf = perf[
            (perf["clf_name"] == META_CLF_NAMES.DNA_COMBINED) & ## Only use DNA combined perfs as RNA combined misses some cancer types
            (perf["cancer_type"].isin(self.predictions.class_columns)) &
            (perf["feat_name"].isin(["n_total", "recall", "precision"]))
        ]

        ## Check if performance stats has all required cancer types
        if perf["cancer_type"].unique().__len__() != self.predictions.class_columns.__len__():
            self.logger.error("`self.cv_performance` is missing some classes found in `self.predictions`")
            raise LookupError

        perf = perf.reindex(columns=self.COLUMN_NAMES)

        perf["data_type"] = "cv_performance"
        perf["clf_group"] = CLF_GROUPS.from_clf_names(perf["clf_name"])

        return perf

    ## Merge ================================
    @staticmethod
    def add_ranks(vis_data: pd.DataFrame) -> None:
        grouper = vis_data.groupby(["sample_id", "data_type", "clf_name", "feat_name"], dropna=False)

        vis_data["rank"] = grouper["data_value"]\
            .rank(method="first", ascending=False, na_option="bottom")\
            .astype(int)\
            .values

        ## Rename rank names to be in ascending order
        rank_groups = grouper.ngroup()
        rank_groups_uniq = rank_groups.unique()
        rank_group_remapping = pd.Series(
            range(rank_groups_uniq.max()+1),
            index=rank_groups.unique()
        )
        rank_groups = rank_groups.replace(rank_group_remapping)

        vis_data["rank_group"] = rank_groups

    def build(self) -> CuppaVisData:

        sample_ids = self.predictions.sample_ids

        if len(sample_ids) == 1:
            if self.verbose:
                self.logger.info("Building vis data for sample '%s'" % sample_ids[0])
        elif self.sample_id is not None:
            if self.verbose:
                self.logger.info("Building vis data for sample '%s' from %i samples" % (self.sample_id, len(sample_ids)))
            self.predictions = self.predictions.get_samples(self.sample_id)
        else:
            if self.verbose:
                self.logger.info("Building vis data for %i samples", len(sample_ids))

        vis_data = pd.concat([
            self.get_probs(),
            self.get_signatures(),
            self.get_top_feat_contrib(),
            self.parse_cv_performance()
        ])

        self.add_ranks(vis_data)
        vis_data.reset_index(drop=True, inplace=True)

        return CuppaVisData.from_data_frame(vis_data)


class CuppaVisData(pd.DataFrame, LoggerMixin):
    def __init__(self, vis_data: pd.DataFrame, *args, **kwargs):
        super().__init__(vis_data, *args, **kwargs)

    @property
    def _constructor(self):
        return CuppaVisData

    def to_tsv(self, path: str, verbose: bool = False):
        if verbose:
            self.logger.info("Writing vis_data to: " + path)

        self.to_csv(path, sep="\t", index=False)

    @staticmethod
    def from_data_frame(df: pd.DataFrame) -> CuppaVisData:
        required_columns = [
            "sample_id", "data_type", "clf_group", "clf_name",
            "feat_name", "feat_value", "cancer_type", "data_value",
            "rank", "rank_group"
        ]
        check_required_columns(df, required_columns=required_columns)

        return CuppaVisData(df)

    @classmethod
    def from_tsv(cls, path: str) -> CuppaVisData:
        df = pd.read_table(path)
        return cls.from_data_frame(df)

    def drop_cv_performance_data(self) -> CuppaVisData:
        return self[self["data_type"] != "cv_performance"]


class CuppaVisPlotter(LoggerMixin):
    def __init__(
        self,
        vis_data: CuppaVisData,
        plot_path: str,
        verbose: bool = True
    ):
        self.vis_data = vis_data
        self.plot_path = os.path.expanduser(plot_path)
        self.verbose = verbose

        self.vis_data_path: str | None = None

        self._check_number_of_samples()
        self._check_plot_path_extension()

    @classmethod
    def from_tsv(cls, path: str, **kwargs) -> CuppaVisPlotter:
        ## This method exists to avoid writing a temporary vis data file if a vis data file already exists
        vis_data = CuppaVisData.from_tsv(path)
        plotter = cls(vis_data=vis_data, **kwargs)
        plotter.vis_data_path = path
        return plotter

    def _check_number_of_samples(self):
        sample_ids = self.vis_data["sample_id"].dropna().unique()

        max_samples = 25
        if len(sample_ids) > max_samples:
            self.logger.error("Plotting predictions for >", max_samples, " is not supported")
            raise RuntimeError

    def _check_plot_path_extension(self):
        if not self.plot_path.endswith((".pdf", ".png")):
            self.logger.error("`plot_path` must end with .pdf or .png")
            raise ValueError

    @property
    def _tmp_vis_data_path(self) -> str:
        return os.path.join(os.path.dirname(self.plot_path), "cuppa.vis_data.tmp.tsv")

    def _write_tmp_vis_data(self):
        if self.verbose:
            self.logger.debug("Writing temporary vis data: " + self._tmp_vis_data_path)
        self.vis_data.to_tsv(self._tmp_vis_data_path)

    def _remove_tmp_vis_data(self):
        if os.path.exists(self._tmp_vis_data_path):
            os.remove(self._tmp_vis_data_path)

            if self.verbose:
                self.logger.debug("Removing temporary vis data: " + self._tmp_vis_data_path)

    def plot(self) -> None:

        try:
            if self.vis_data_path is None or not os.path.exists(self.vis_data_path):
                self._write_tmp_vis_data()
                self.vis_data_path = self._tmp_vis_data_path

            executor = RscriptExecutor([
                RSCRIPT_PLOT_PREDICTIONS_PATH,
                self.vis_data_path,
                self.plot_path
            ])
            executor.run()

        finally:
            self._remove_tmp_vis_data()
