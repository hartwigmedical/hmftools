from __future__ import annotations

import os
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
import tempfile

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
        verbose: bool = True
    ):
        self.predictions = predictions
        self.sample_id = sample_id
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

    def _subset_feat_contrib_by_feat_pattern(self, pattern: str) -> CuppaPrediction:
        return self.feat_contrib[
            self.feat_contrib.index.get_level_values("feat_name").str.contains(pattern)
        ]

    def _get_top_drivers(
        self,
        min_driver_contrib: Optional[float] = 0.2,
        min_driver_likelihood: Optional[float] = 0.2,
    ) -> CuppaPrediction:
        feat_contrib = self._subset_feat_contrib_by_feat_pattern("driver")

        index = feat_contrib.index.to_frame(index=False)
        row_maxes = (feat_contrib.max(axis=1) >= min_driver_contrib).values
        selected_rows = (row_maxes >= min_driver_contrib) | (index["feat_value"] >= min_driver_likelihood)

        return feat_contrib[selected_rows.values]

    def _get_existing_features(self, feat_pattern: str):
        feat_contrib = self._subset_feat_contrib_by_feat_pattern(feat_pattern)
        max_contrib_by_row = feat_contrib.max(axis=1)
        feat_values = feat_contrib.index.get_level_values("feat_value")
        return feat_contrib[(max_contrib_by_row > 0) & (feat_values > 0)]

    def _get_existing_fusions(self) -> CuppaPrediction:
        return self._get_existing_features("fusion")

    def _get_existing_viruses(self) -> CuppaPrediction:
        return self._get_existing_features("virus")

    def get_top_feat_contrib(
        self,
        min_driver_contrib: Optional[float] = 0.2,
        min_driver_likelihood: Optional[float] = 0.2,
    ) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting top event features")

        wide = self.predictions.__class__.concat([
            self._subset_feat_contrib_by_feat_pattern("trait"),
            self._subset_feat_contrib_by_feat_pattern("tmb"),

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

        wide = self.predictions.get_data_types("cv_performance")
        wide = wide[
            wide.index.get_level_values("feat_name").isin(required_metrics) &
            (wide.index.get_level_values("clf_name") == clf_name)
        ]

        ## Force metric order
        long = wide.wide_to_long()
        long["feat_name"] = pd.Categorical(long["feat_name"], required_metrics)
        long = long.sort_values("feat_name")

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
    def __init__(
        self,
        vis_data: CuppaVisData,
        plot_path: str,
        vis_data_path: Optional[str] = None,
        verbose: bool = True
    ):
        self.vis_data = vis_data
        self.plot_path = os.path.expanduser(plot_path)
        self._check_plot_path()

        self.verbose = verbose

        if vis_data_path is None:
            self._using_tmp_vis_data_path = True
            self.vis_data_path = self._get_tmp_vis_data_path()
        else:
            self.vis_data_path = vis_data_path
            self._using_tmp_vis_data_path = False

    def _check_plot_path(self) -> None:

        if not self.plot_path.endswith((".pdf", ".png", ".jpg")):
            self.logger.error("`plot_path` must end with .pdf, .png, or .jpg")
            raise ValueError

    def _get_tmp_vis_data_path(self) -> str:
        ## If no prediction_path is specified, write to tmp location so that it can be passed to R
        vis_data_path = os.path.join(tempfile.gettempdir(), "vis_data.tsv")

        if self.verbose:
            self.logger.info("`vis_data_path` not specified. Using temporary path: " + vis_data_path)

        return vis_data_path

    def remove_tmp_vis_data(self) -> None:
        if not self._using_tmp_vis_data_path:
            return None

        if self.verbose:
            self.logger.info("Removing temporary vis data")

    def run_r_script(self) -> None:
        command = f"Rscript --vanilla {RSCRIPT_PLOT_PREDICTIONS_PATH} {self.vis_data_path} {self.plot_path}"

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
        if not os.path.exists(self.vis_data_path):
            self.vis_data.to_tsv(self.vis_data_path)
        else:
            self.logger.info("Using existing vis_data file at: " + self.vis_data_path)

        self.run_r_script()
        self.remove_tmp_vis_data()
