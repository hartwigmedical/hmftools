from __future__ import annotations

import os
from functools import cached_property
from pprint import pformat
from typing import Optional

import pandas as pd
from sklearn.model_selection import StratifiedKFold

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.runners.args import DEFAULT_RUNNER_ARGS
from cuppa.compose.pipeline import PipelineCrossValidator
from cuppa.logger import LoggerMixin, reset_logging_basic_config
from cuppa.performance.confusion_matrix import ConfusionMatrix
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, FeatureLoader, CuppaFeatures
from cuppa.sample_data.sample_metadata import SampleMetadata, TrainingSampleSelector
from cuppa.performance.performance_stats import PerformanceStats


class TrainingRunner(LoggerMixin):

    def __init__(
        self,

        features_path: str,
        output_dir: str,
        metadata_path: str,

        using_old_features_format: bool = DEFAULT_RUNNER_ARGS.using_old_features_format,

        genome_version: int = DEFAULT_RUNNER_ARGS.genome_version,
        min_samples_with_rna: int = DEFAULT_RUNNER_ARGS.min_samples_with_rna,
        excl_classes: str | list[str] = DEFAULT_RUNNER_ARGS.excl_classes,
        excl_chroms: str | list[str] = DEFAULT_RUNNER_ARGS.excl_chroms,

        fusion_overrides_path=DEFAULT_RUNNER_ARGS.fusion_overrides_path,

        skip_cv: bool = DEFAULT_RUNNER_ARGS.skip_cv,
        cv_folds: int = DEFAULT_RUNNER_ARGS.cv_folds,
        cache_training: bool = DEFAULT_RUNNER_ARGS.cache_training,
        n_jobs: int = DEFAULT_RUNNER_ARGS.n_jobs,

        log_to_file: bool = DEFAULT_RUNNER_ARGS.log_to_file,
        log_path: Optional[str] = DEFAULT_RUNNER_ARGS.log_path
    ):
        ## Args --------------------------------
        ## Data
        self.features_path = features_path
        self.output_dir = output_dir
        self.metadata_path = metadata_path

        ## Settings
        self.using_old_features_format = using_old_features_format

        self.genome_version = genome_version
        self.excl_chroms = excl_chroms
        self.excl_classes = excl_classes
        self.min_samples_with_rna = min_samples_with_rna

        self.fusion_overrides_path = fusion_overrides_path

        ## CV
        self.skip_cv = skip_cv
        self.cv_folds = cv_folds
        self.cache_training = cache_training
        self.n_jobs = n_jobs

        ## Logging
        self.log_to_file = log_to_file
        self.log_path = log_path
        self.set_log_path()

        ## Attributes assigned at run time --------------------------------
        self.sample_metadata: SampleMetadata = None
        self.training_samples: dict[str, pd.Index] = None
        self.X: CuppaFeatures = None
        self.y: pd.Series = None
        self.y_split: pd.Series = None

        self.cross_validator: PipelineCrossValidator = None
        self.cv_predictions: CuppaPrediction = None
        self.cv_pred_summ: CuppaPredSummary = None
        self.cv_performance: PerformanceStats = None
        self.cv_performance_by_prob_bin: PerformanceStats = None

        self.cuppa_classifier: CuppaClassifier = None


    def __repr__(self) -> str:
        return self.__class__.__name__ + " args:\n" + pformat(vars(self))

    def set_log_path(self) -> None:
        if not self.log_to_file:
            return None

        log_path = self.log_path
        if log_path is None:
            log_path = os.path.join(self.output_dir, "run.log")

        reset_logging_basic_config(filename=log_path)

    ## Load data ================================
    def load_sample_metadata(self) -> None:
        ## TODO: update to from_tsv() when the new round of training comes
        self.sample_metadata = SampleMetadata.from_csv(self.metadata_path)

    def get_training_samples(self) -> None:

        self.sample_selector = TrainingSampleSelector(
            metadata=self.sample_metadata,
            excl_classes=self.excl_classes,
            min_samples_with_rna=self.min_samples_with_rna,
            verbose=True,
        )

        self.training_samples = self.sample_selector.selected_samples

    def write_training_samples(self) -> None:
        metadata_out_path = os.path.join(self.output_dir, "metadata.tsv.gz")
        self.logger.info("Writing training metadata to: " + metadata_out_path)
        self.sample_selector.metadata.to_tsv(metadata_out_path)

    def get_y(self) -> None:
        self.y = self.sample_metadata.loc[self.training_samples["dna"], "CancerSubtype"]

    def _get_y_split(self) -> None:
        ## Get label for cross-validation splits based on both cancer type and rna read length
        rna_read_length = self.sample_metadata.loc[self.y.index, "RnaReadLength"]
        y_split = self.y + "__" + rna_read_length.astype(str)
        self.y_split = y_split

    def get_X(self) -> None:

        if not self.using_old_features_format:
            loader = FeatureLoader(self.features_path)
            X = loader.load()
            self.X = X
            return None
            ## TODO: add DNA and RNA sample selection


        paths = CuppaFeaturesPaths.from_dir(self.features_path, file_format="old")

        loader = FeatureLoaderOld(
            paths=paths,
            genome_version=self.genome_version,
            excl_chroms=self.excl_chroms,
            verbose=True
        )

        X_dna = loader.load_dna_features()
        X_rna = loader.load_rna_features()

        X_dna = X_dna.loc[self.training_samples["dna"]]
        X_rna = X_rna.loc[self.training_samples["rna"]]

        X = pd.concat([X_dna, X_rna], axis=1)

        self.X = X

    def load_data(self) -> None:
        self.load_sample_metadata()
        self.get_training_samples()

        self.get_y()
        self._get_y_split()

        self.get_X()

    ## CV ================================
    def cv_fit(self) -> None:

        self._cv_dir = None
        if self.cache_training:
            self._cv_dir = os.path.join(self.output_dir, "cv")

        self.cross_validator = PipelineCrossValidator(
            pipeline=CuppaClassifier(fusion_overrides_path=self.fusion_overrides_path),
            X=self.X,
            y=self.y,
            y_split=self.y_split,

            out_dir=self._cv_dir,

            cv=StratifiedKFold(n_splits=self.cv_folds, random_state=0, shuffle=True),
            n_jobs=self.n_jobs
        )

        self.cross_validator.fit(cache_training=self.cache_training)

    def get_cv_predictions(self) -> None:
        self.cv_predictions = self.cross_validator.apply_on_test_sets("predict", n_jobs=1)

    def get_cv_pred_summ(self) -> None:
        predictions = self.cv_predictions
        self.cv_pred_summ = predictions.summarize(actual_classes=self.y, show_extra_info=True, verbose=True)

    def get_cv_performance(self) -> None:
        self.cv_performance = self.cv_pred_summ.performance()

    def get_cv_performance_by_prob_bin(self) -> None:
        self.cv_performance_by_prob_bin = self.cv_pred_summ.performance_by_prob_bin()

    def plot_confusion_matrices(self) -> None:
        confusion_dir = os.path.join(self.cv_report_dir, "confusion")
        os.makedirs(confusion_dir, exist_ok=True)

        for clf_name in self.cv_pred_summ.clf_names:
            pred_summ_clf = self.cv_pred_summ.query(f"clf_name=='{clf_name}'")
            confusion_matrix = ConfusionMatrix(pred_summ_clf, clf_name=clf_name)
            confusion_matrix.plot(os.path.join(confusion_dir, f"confusion.{clf_name}.pdf"))

    @cached_property
    def cv_report_dir(self) -> str:
        path = os.path.join(self.output_dir, "cv/report")
        os.makedirs(path, exist_ok=True)
        return path

    def make_cv_report(self) -> None:

        self.get_cv_predictions()
        self.get_cv_pred_summ()
        self.get_cv_performance()
        self.get_cv_performance_by_prob_bin()

        self.cv_pred_summ.to_tsv(os.path.join(self.cv_report_dir, "pred_summ.tsv.gz"), verbose=True)

        self.cv_performance.to_tsv(os.path.join(self.cv_report_dir, "perf.tsv"), verbose=True)
        self.cv_performance_by_prob_bin.to_tsv(os.path.join(self.cv_report_dir, "perf_by_prob_bin.tsv"), verbose=True)

        self.plot_confusion_matrices()

        self.cv_predictions.to_tsv(os.path.join(self.cv_report_dir, "predictions.tsv.gz"), verbose=True)

    ## Final model ================================
    def train_final_model(self) -> None:

        fit_cache_dir = None
        if self.cache_training:
            fit_cache_dir = os.path.join(self.output_dir, "model_cache")

        self.cuppa_classifier = CuppaClassifier(
            fit_cache_dir=fit_cache_dir,
            fusion_overrides_path = self.fusion_overrides_path,
            n_jobs=self.n_jobs
        )

        self.cuppa_classifier.fit(self.X, self.y)

        if self.cv_performance is not None:
            self.cuppa_classifier.cv_performance = self.cv_performance
        else:
            self.logger.warning("Could not add CV performance stats to CuppaClassifier")

    def export_final_model(self) -> None:
        path = os.path.join(self.output_dir, "cuppa_classifier.pickle.gz")
        self.cuppa_classifier.to_file(path)

    ## Run ================================
    def run(self) -> None:

        self.load_data()
        self.write_training_samples()

        if not self.skip_cv:
            self.logger.info("Performing cross-validation")
            self.cv_fit()

            self.logger.info("Making cross-validation report")
            self.make_cv_report()

        self.logger.info("Training final model")
        self.train_final_model()
        self.export_final_model()