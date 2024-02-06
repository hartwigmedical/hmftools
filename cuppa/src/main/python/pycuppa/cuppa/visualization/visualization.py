from __future__ import annotations

import os
from subprocess import Popen, PIPE, STDOUT, CalledProcessError

import numpy as np
import pandas as pd
from functools import cached_property

from cuppa.constants import RSCRIPT_PLOT_PREDICTIONS_PATH
from cuppa.logger import LoggerMixin
from cuppa.misc.utils import check_required_columns

from typing import TYPE_CHECKING, Optional


if TYPE_CHECKING:
    from cuppa.classifier.cuppa_prediction import CuppaPrediction


class CuppaVisDataBuilder(LoggerMixin):
    def __init__(
        self,
        predictions: CuppaPrediction,
        sample_id: Optional[str | int] = None,
        require_all_feat_types: bool = True,
        verbose: bool = True
    ):
        self.predictions = predictions
        self.sample_id = sample_id
        self.require_all_feat_types = require_all_feat_types
        self.verbose = verbose

        self._get_predictions_one_sample()

    def _get_predictions_one_sample(self) -> CuppaPrediction:
        if not self.predictions.is_multi_sample:
            return self.predictions

        if self.sample_id is None:
            self.logger.error("Cannot build vis data of many samples. Please provide `sample_id` when `predictions` contains multiple samples")
            raise Exception

        predictions = self.predictions.get_samples(self.sample_id)

        if self.predictions.has_cv_performance:
            predictions = predictions.__class__.concat([
                predictions,
                self.predictions.get_data_types("cv_performance")
            ])
        else:
            self.logger.warning("`cv_performance` data missing from `predictions`")

        self.predictions = predictions

    def _check_snv_counts_feature_exists(self):
        FEAT_NAME_SNV_COUNT = "event.tmb.snv_count"
        if FEAT_NAME_SNV_COUNT not in self.predictions.index.get_level_values("feat_name"):
            self.logger.error(
                "Feature contributions / values for '%s' do not exist. "
                "This is required for calculating relative mutational signature counts"
            )
            raise LookupError

    ## Probs ================================
    def get_probs(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.info("Getting probabilities")

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
                self.logger.error("No features found with pattern '%s'" % pattern)
                raise LookupError
            else:
                self.logger.warning("No features found with pattern '%s'. Returning None" % pattern)
                return None

        return feat_contrib

    def _get_top_drivers(
        self,
        min_driver_contrib: Optional[float] = 0.2,
        min_driver_likelihood: Optional[float] = 0.2,
    ) -> CuppaPrediction | None:
        feat_contrib = self._subset_feat_contrib_by_feat_pattern("driver")

        if feat_contrib is None:
            return None

        index = feat_contrib.index.to_frame(index=False)
        row_maxes = (feat_contrib.max(axis=1) >= min_driver_contrib).values
        selected_rows = (row_maxes >= min_driver_contrib) | (index["feat_value"] >= min_driver_likelihood)

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

    def get_top_feat_contrib(
        self,
        min_driver_contrib: Optional[float] = 0.2,
        min_driver_likelihood: Optional[float] = 0.2,
    ) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting top event features")

        wide = self.predictions.__class__.concat([  ## Indirect call to CuppaPrediction.concat() to avoid circular import
            self._subset_feat_contrib_by_feat_pattern("tmb"),
            self._subset_feat_contrib_by_feat_pattern("trait"),

            self._get_existing_fusions(),
            self._get_existing_viruses(),

            self._get_top_drivers(
                min_driver_contrib=min_driver_contrib,
                min_driver_likelihood=min_driver_likelihood
            ),

            self._subset_feat_contrib_by_feat_pattern("sv")
        ])

        long = wide.wide_to_long()

        ## Calculate odds
        long["data_value"] = np.exp(1) ** long["data_value"]

        ## Strip clf_name from feat_name
        long["feat_name"] = long["feat_name"].str.split(".", n=1).map(lambda x: x[1])

        return long

    ## CV performance ================================
    def _get_cv_performance_one_clf(self, clf_name: str = "dna_combined") -> pd.DataFrame | None:

        if not self.predictions.has_cv_performance:
            return None

        if self.verbose:
            self.logger.debug("Getting cross-validation performance")

        required_metrics = ["n_total", "recall", "precision"]
        cancer_type_order = self.predictions.class_metadata.index

        wide = self.predictions.get_data_types("cv_performance")
        wide = wide[
            wide.index.get_level_values("feat_name").isin(required_metrics) &
            (wide.index.get_level_values("clf_name") == clf_name)
        ]

        ## Force metric and cancer type order
        long = wide.wide_to_long()
        long["feat_name"] = pd.Categorical(long["feat_name"], required_metrics)
        long["cancer_type"] = pd.Categorical(long["cancer_type"], cancer_type_order)
        long = long.sort_values(["feat_name","cancer_type"])

        ## Rename
        long["feat_name"] = long["feat_name"].replace(dict(n_total="n_samples"))

        return long

    def get_cv_performance(self) -> pd.DataFrame | None:
        ## Only use DNA combined perfs as RNA combined misses some cancer types
        return self._get_cv_performance_one_clf(clf_name="dna_combined")

    ## Merge ================================
    @staticmethod
    def _add_rank(vis_data: pd.DataFrame) -> None:
        grouper = vis_data.groupby(["data_type", "clf_name", "feat_name"], dropna=False)

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

        sample_id = self.predictions.sample_ids[0]
        if self.verbose:
            self.logger.info("Building vis data for sample: " + sample_id)

        vis_data = pd.concat([
            self.get_probs(),
            self.get_signatures(),
            self.get_top_feat_contrib(),
            self.get_cv_performance()
        ])

        self._add_rank(vis_data)

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
    def __init__(self, vis_data: CuppaVisData, plot_path: str, verbose: bool = True):
        self.vis_data = vis_data
        self._check_vis_data_has_one_sample()

        self.plot_path = os.path.expanduser(plot_path)
        self._check_plot_path_extension()

        self.verbose = verbose

    def _check_vis_data_has_one_sample(self):
        sample_ids = self.vis_data["sample_id"].dropna().unique()
        if len(sample_ids) > 1:
            self.logger.error("Cannot plot visualization when `vis_data` contains multiple samples")
            raise Exception

    def _check_plot_path_extension(self):
        if not self.plot_path.endswith((".pdf", ".png", ".jpg")):
            self.logger.error("`plot_path` must end with .pdf, .png, or .jpg")
            raise ValueError

    @property
    def tmp_vis_data_path(self):
        return os.path.join(os.path.dirname(self.plot_path), "cuppa.vis_data.tmp.tsv")

    def write_tmp_vis_data(self):
        if self.verbose:
            self.logger.debug("Writing vis data to temporary path: " + self.tmp_vis_data_path)
        self.vis_data.to_tsv(self.tmp_vis_data_path)

    def remove_tmp_vis_data(self):
        if self.verbose:
            self.logger.debug("Removing temporary vis data at: " + self.tmp_vis_data_path)
        os.remove(self.tmp_vis_data_path)

    def plot_in_r(self, vis_data_path: str, plot_path: str) -> None:
        command = f"Rscript --vanilla {RSCRIPT_PLOT_PREDICTIONS_PATH} {vis_data_path} {plot_path}"

        if self.verbose:
            self.logger.info("Running shell command: '%s'" % command)

        with Popen(command, shell=True, stdout=PIPE, stderr=STDOUT) as process:
            for line in iter(process.stdout.readline, b''):
                r_stderr = line.decode("utf-8").strip()
                self.logger.error("[R process] " + r_stderr)

        return_code = process.poll()
        if return_code:
            raise CalledProcessError(return_code, command)

    def plot(self) -> None:
        self.write_tmp_vis_data()

        try:
            self.plot_in_r(vis_data_path=self.tmp_vis_data_path, plot_path=self.plot_path)
            self.remove_tmp_vis_data()
        except CalledProcessError:
            self.remove_tmp_vis_data()

    @classmethod
    def plot_from_tsv(cls, vis_data_path: str, plot_path: str, verbose: bool = True):
        ## This utility method exists to avoid writing a temporary vis data file if a vis data file already exists
        vis_data = CuppaVisData.from_tsv(vis_data_path)
        plotter = cls(vis_data=vis_data, plot_path=plot_path, verbose=verbose)
        plotter.plot()

