from __future__ import annotations

import os
from functools import cached_property
from pprint import pformat
from typing import Optional

import pandas as pd

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.runners.args import DEFAULT_RUNNER_ARGS
from cuppa.logger import LoggerMixin, initialize_logging
from cuppa.sample_data.cuppa_features import CuppaFeaturesLoader
from cuppa.visualization.visualization import CuppaVisData, CuppaVisPlotter, CuppaVisDataBuilder


class PredictionRunner(LoggerMixin):
    def __init__(
        self,
        classifier_path: str,
        features_path: str,
        output_dir: str,
        sample_id: str | None = None,
        compress_tsv_files: bool = False,
        cv_predictions_path: str = None,
        cv_predictions: CuppaPrediction | None = None,
        clf_group: str = DEFAULT_RUNNER_ARGS.clf_group,
        log_to_file: bool = DEFAULT_RUNNER_ARGS.log_to_file,
        log_path: Optional[str] = DEFAULT_RUNNER_ARGS.log_path,
        log_format: str = DEFAULT_RUNNER_ARGS.log_format
    ):
        ## Args --------------------------------
        self.features_path = features_path
        self.output_dir = output_dir
        self.sample_id = sample_id
        self.compress_tsv_files = compress_tsv_files
        self.classifier_path = classifier_path

        self.cv_predictions_path = cv_predictions_path
        self.cv_predictions = cv_predictions
        self.clf_group = clf_group

        self.log_to_file = log_to_file
        self.log_path = log_path
        self.log_format = log_format
        self.set_up_logging()

        ## Attributes assigned at run time --------------------------------
        self.X: pd.DataFrame = None
        self.predictions: CuppaPrediction = None
        self.pred_summ: CuppaPredSummary = None
        self.vis_data: CuppaVisData = None

    def __repr__(self) -> str:
        return self.__class__.__name__ + " args:\n" + pformat(vars(self))

    def set_up_logging(self) -> None:
        if self.log_to_file and self.log_path is None:
            self.log_path = os.path.join(self.output_dir, "predict.log")

        initialize_logging(filename=self.log_path, format=self.log_format)

    @cached_property
    def cuppa_classifier(self) -> CuppaClassifier:
        return CuppaClassifier.from_file(self.classifier_path)

    def get_X(self) -> None:

        loader = CuppaFeaturesLoader(self.features_path, sample_id=self.sample_id)
        X = loader.load()

        X = self.cuppa_classifier.fill_missing_cols(X)

        self.X = X

    @property
    def n_samples(self) -> int:
        return self.X.shape[0]

    def get_predictions(self) -> None:

        if self.cv_predictions_path is not None:
            self.logger.info("Loading cross-validation predictions from: " + self.cv_predictions_path)
            self.cv_predictions = CuppaPrediction.from_tsv(self.cv_predictions_path)

        if self.cv_predictions is None:
            self.logger.info("Computing predictions for %i samples" % self.n_samples)
            self.predictions = self.cuppa_classifier.predict(self.X)
            return

        is_cv_sample = self.X.index.isin(self.cv_predictions.sample_ids)
        sample_ids_new = self.X.index[~is_cv_sample]
        sample_ids_cv = self.X.index[is_cv_sample]

        predictions = []

        if len(sample_ids_new) > 0:
            self.logger.info("Computing predictions for %i / %i samples" % (len(sample_ids_new), self.n_samples))
            X_new = self.X.loc[sample_ids_new]
            predictions.append(self.cuppa_classifier.predict(X_new))

        if len(sample_ids_cv) > 0:
            predictions.append(self.cv_predictions.loc[sample_ids_cv])
            self.logger.info("Using pre-computed cross-validation predictions for %i / %i samples" % (len(sample_ids_cv), self.n_samples))

        predictions = CuppaPrediction.concat(predictions)
        predictions = predictions.loc[self.X.index]

        if self.clf_group.lower() == "all":
            pass
        elif self.clf_group.lower() == "dna":
            predictions = predictions.subset_probs_by_clf_groups("dna")
        else:
            self.logger.error("`clf_group` must be 'all' or 'dna'")
            raise ValueError

        self.predictions = predictions

    def get_pred_summ(self) -> None:
        self.pred_summ = self.predictions.summarize(show_extra_info=True, verbose=True)

    def get_vis_data(self) -> None:
        builder = CuppaVisDataBuilder(
            predictions = self.predictions,
            sample_id = self.sample_id,
            cv_performance = self.cuppa_classifier.cv_performance
        )

        self.vis_data = builder.build()

    def add_filename_affixes(self, filename: str):
        if self.sample_id is not None:
            filename = self.sample_id + "." + filename

        if self.compress_tsv_files and filename.endswith(".tsv"):
            filename += ".gz"

        return filename

    @property
    def vis_data_path(self) -> str:
        filename = self.add_filename_affixes("cuppa.vis_data.tsv")
        return os.path.join(self.output_dir, filename)

    @property
    def plot_path(self) -> str:

        if self.n_samples == 1 or self.sample_id is not None:
            ## For single sample, use png so that it can be inserted into the ORANGE pdf
            filename = "cuppa.vis.png"
        else:
            ## For multi sample, use pdf so that the viz for all samples can be contained in one pdf
            filename = "cuppa.vis.pdf"

        filename = self.add_filename_affixes(filename)

        return os.path.join(self.output_dir, filename)

    @property
    def pred_summ_path(self) -> str:
        filename = self.add_filename_affixes("cuppa.pred_summ.tsv")
        return os.path.join(self.output_dir, filename)

    def run(self) -> None:
        self.get_X()

        if self.n_samples == 1 and self.sample_id is None:
            self.logger.info("`sample_id` must be provided when predicting on one sample")
            raise ValueError

        self.get_predictions()
        self.get_pred_summ()
        self.get_vis_data()

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            self.logger.info("Created output dir: " + self.output_dir)

        self.pred_summ.to_tsv(self.pred_summ_path, verbose=True)
        self.vis_data.to_tsv(self.vis_data_path, verbose=True)

        CuppaVisPlotter.from_tsv(path=self.vis_data_path, plot_path=self.plot_path, verbose=True).plot()
