import argparse
import os
import sys
from functools import cached_property
from pprint import pformat
from typing import Any

import pandas as pd
from sklearn.model_selection import StratifiedKFold

from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary
from cuppa.compose.pipeline import PipelineCrossValidator
from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH, DEFAULT_CUPPA_CLASSIFIER_PATH
from cuppa.logger import LoggerMixin, reset_logging_basic_config
from cuppa.performance.confusion_matrix import ConfusionMatrix
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, FeatureLoaderNew, CuppaFeatures
from cuppa.sample_data.sample_metadata import SampleMetadata, TrainingSampleSelector
from cuppa.performance.performance_stats import PerformanceStats
from cuppa.visualization.visualization import CuppaVisData, CuppaVisPlotter


class DEFAULT_RUNNER_ARGS:
    classifier_path: str = DEFAULT_CUPPA_CLASSIFIER_PATH

    output_prefix: str = None
    compress_tsv_files: bool = False
    using_old_features_format: bool = False

    genome_version: int = 37
    min_samples_with_rna: int = 5
    excl_chroms: str | list[str] = ["ChrY", "Y"]
    excl_classes: str | list[str] = ["_Other","_Unknown"]

    fusion_overrides_path = DEFAULT_FUSION_OVERRIDES_PATH

    skip_cv: bool = False
    cv_folds: int = 10
    cache_training: bool = True
    n_jobs: int = 1

    log_to_file: bool = False
    log_path: str | None = None


class RunnerArgParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser()

    @staticmethod
    def comma_sep_str_to_list(string) -> list[str]:
        return [s.strip() for s in string.split(",")]

    ## Arguments ================================
    ## I/O --------------------------------
    def add_features_path(self) -> None:
        self.parser.add_argument(
            "--features_path",
            type=str, required=True,
            help="[Old features format] Directory with the input features csv files. [New features format] Path to the input features tsv file"
        )

    def add_output_dir(self) -> None:
        self.parser.add_argument(
            "--output_dir",
            type=str, required=True, help="Output directory"
        )

    def add_metadata_path(self) -> None:
        self.parser.add_argument(
            "--metadata_path",
            type=str, required=True, help="Path to the metadata file with cancer type labels per sample"
        )

    def add_classifier_path(self) -> None:
        self.parser.add_argument(
            "--classifier_path",
            help="Path to the pickle file of the classifier."
                 "If not provided, the default classifier path will be used (pycuppa/resources/cuppa_classifier.pickle.gz)"
        )

    def add_sample_id(self) -> None:
        self.parser.add_argument(
            "--sample_id",
            help="If provided, will prepend `sample_id` to the output filenames"
        )

    def add_compress_tsv_files(self):
        self.parser.add_argument(
            "--compress_tsv_files", dest="compress_tsv_files", action="store_true",
            help="Compress tsv files with gzip? (will add .gz to the file extension)"
        )

    def add_using_old_features_format(self) -> None:
        self.parser.add_argument("--using_old_features_format", action="store_true", help="Use old input features format?")

    ## Features/ metadata options --------------------------------
    def add_genome_version(self) -> None:
        self.parser.add_argument(
            "--genome_version", type=int, default=DEFAULT_RUNNER_ARGS.genome_version, choices=[37, 38],
            help="Genome version. Can be 37 or 38. Used for getting the gen_pos bin names. Default=" + str(DEFAULT_RUNNER_ARGS.genome_version)
        )

    def add_excl_chroms(self) -> None:
        self.parser.add_argument(
            "--excl_chroms", type=self.comma_sep_str_to_list,
            default=DEFAULT_RUNNER_ARGS.excl_chroms,
            help="Comma separated list of chromosomes to exclude from the gen_pos matrix."
                 "E.g. 'chrY' or 'chr1,chr2'. Default: "
                 ",".join(DEFAULT_RUNNER_ARGS.excl_chroms)
        )

    def add_min_samples_with_rna(self) -> None:
        self.parser.add_argument(
            "--min_samples_with_rna", type=int, default=DEFAULT_RUNNER_ARGS.min_samples_with_rna,
            help="Minimum number of samples with RNA in each cancer subtype. "
                 "If the cancer subtype has fewer samples with RNA than this value, the cancer subtype will be excluded "
                 "from training. Default: " + str(DEFAULT_RUNNER_ARGS.min_samples_with_rna)
        )

    def add_excl_classes(self) -> None:
        self.parser.add_argument(
            "--excl_classes", type=self.comma_sep_str_to_list,
            default=DEFAULT_RUNNER_ARGS.excl_classes,
            help="Comma separated list of cancer subtypes to exclude from training. "
                 "E.g. 'Breast' or 'Breast,Lung'. Default: "
                 ",".join(DEFAULT_RUNNER_ARGS.excl_classes)
        )

    def add_fusion_overrides_path(self) -> None:
        self.parser.add_argument(
            "--fusion_overrides_path", type=str, default=DEFAULT_RUNNER_ARGS.fusion_overrides_path,
            help="Path to the fusion overrides tsv file"
        )

    ## Cross-validation / training --------------------------------
    def add_cv_folds(self) -> None:
        self.parser.add_argument(
            "--cv_folds", type=int, default=DEFAULT_RUNNER_ARGS.cv_folds,
            help="Number of cross-validation folds. Default: " + str(DEFAULT_RUNNER_ARGS.cv_folds)
        )

    def add_skip_cv(self) -> None:
        self.parser.add_argument("--skip_cv", action="store_true", help="Skip cross-validation?")

    def add_no_cache_training(self):
        self.parser.add_argument(
            "--no_cache_training", dest="cache_training", action="store_false",
            help="Don't cache models during cross-validation and training of final model?"
        )

    def add_n_jobs(self) -> None:
        self.parser.add_argument(
            "--n_jobs", type=int, default=DEFAULT_RUNNER_ARGS.n_jobs,
            help="Number of threads to use during cross-validation (parallelized cross-validation folds) "
                 "and final model training (ColumnTransformer -> parallelized sub-/meta-classifier training). "
                 "If -1, all threads are used. Default" + str(DEFAULT_RUNNER_ARGS.n_jobs)
        )

    ## Logging --------------------------------
    def add_log_to_file(self):
        self.parser.add_argument("--log_to_file", dest="log_to_file", action="store_true", help="Output logs to a file?")

    def add_log_path(self) -> None:
        self.parser.add_argument(
            "--log_path", type=str, default=DEFAULT_RUNNER_ARGS.log_path,
            help="Path to log file. Default: $(output_dir)/run.log"
        )

    ## Arg groups ================================
    def get_kwargs_predict(self) -> dict[str, Any]:
        self.parser.description = "PredictionRunner"

        self.add_features_path()
        self.add_output_dir()
        self.add_sample_id()
        self.add_compress_tsv_files()

        self.add_classifier_path()

        self.add_using_old_features_format()
        self.add_genome_version()
        self.add_excl_chroms()
        self.add_log_to_file()
        self.add_log_path()

        args = self.parser.parse_args()
        return args.__dict__

    def get_kwargs_train(self) -> dict[str, Any]:
        self.parser.description = "TrainingRunner"

        self.add_features_path()
        self.add_output_dir()
        self.add_metadata_path()

        self.add_using_old_features_format()
        self.add_genome_version()
        self.add_min_samples_with_rna()
        self.add_excl_classes()
        self.add_excl_chroms()

        self.add_fusion_overrides_path()

        self.add_skip_cv()
        self.add_cv_folds()
        self.add_no_cache_training()
        self.add_n_jobs()

        self.add_log_to_file()
        self.add_log_path()

        args = self.parser.parse_args()
        return args.__dict__


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
        log_path: str | None = DEFAULT_RUNNER_ARGS.log_path
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
            loader = FeatureLoaderNew(self.features_path)
            X = loader.load()
            self.X = X
            return None
            ## TODO: add DNA and RNA sample selection


        paths = CuppaFeaturesPaths.from_dir(self.features_path, basenames_mode="old")

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

    def load_data(self):
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
        self.cv_predictions = self.cross_validator.apply_on_test_sets("predict", probs_only=False, n_jobs=1)  ## TODO: add probs_only as a global argument

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
            confusion_matrix.plot(os.path.join(confusion_dir, f"confusion.{clf_name}.pdf"), verbose=True)

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
            self.logger.warning("Could add CV performance stats to CuppaClassifier")

    def export_final_model(self) -> None:
        path = os.path.join(self.output_dir, "cuppa_classifier.pickle.gz")
        self.cuppa_classifier.to_file(path, verbose=True)

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


class PredictionRunner(LoggerMixin):
    def __init__(
        self,
        features_path: str,
        output_dir: str,
        sample_id: str | None = None,
        compress_tsv_files: bool = False,
        classifier_path: str | None = None,
        using_old_features_format: bool = DEFAULT_RUNNER_ARGS.using_old_features_format,
        genome_version: int = DEFAULT_RUNNER_ARGS.genome_version,
        excl_chroms: str | list[str] = DEFAULT_RUNNER_ARGS.excl_chroms,
        log_to_file: bool = DEFAULT_RUNNER_ARGS.log_to_file,
        log_path: str | None = DEFAULT_RUNNER_ARGS.log_path
    ):
        ## Args --------------------------------
        self.features_path = features_path
        self.output_dir = output_dir
        self.sample_id = sample_id
        self.compress_tsv_files = compress_tsv_files
        self.classifier_path = DEFAULT_RUNNER_ARGS.classifier_path if classifier_path is None else classifier_path

        self.using_old_features_format = using_old_features_format

        self.genome_version = genome_version
        self.excl_chroms = excl_chroms

        self.log_to_file = log_to_file
        self.log_path = log_path
        self.set_log_path()

        ## Attributes assigned at run time --------------------------------
        self.predictions: CuppaPrediction = None
        self.pred_summ: CuppaPredSummary = None
        self.vis_data: CuppaVisData = None

    def __repr__(self) -> str:
        return self.__class__.__name__ + " args:\n" + pformat(vars(self))

    def set_log_path(self) -> None:
        if not self.log_to_file:
            return None

        log_path = self.log_path
        if log_path is None:
            log_path = os.path.join(self.output_dir, "run.log")

        reset_logging_basic_config(filename=log_path)

    @cached_property
    def cuppa_classifier(self) -> CuppaClassifier:

        if self.classifier_path == DEFAULT_RUNNER_ARGS.classifier_path:
            msg = "Loading default classifier from pycuppa resources: " + str(self.classifier_path)
        else:
            msg = "Loading classifier from: " + str(self.classifier_path)

        self.logger.info(msg)

        return CuppaClassifier.from_file(self.classifier_path)

    def get_X(self):

        if self.using_old_features_format:
            paths = CuppaFeaturesPaths.from_dir(self.features_path, basenames_mode="old")
            loader = FeatureLoaderOld(
                paths=paths,
                genome_version=self.genome_version,
                excl_chroms=self.excl_chroms,
                verbose=True
            )

            X = loader.load_features()
        else:
            loader = FeatureLoaderNew(self.features_path)
            X = loader.load()

        X = self.cuppa_classifier.fill_missing_cols(X)

        self.X = X

    def get_predictions(self) -> None:
        self.predictions = self.cuppa_classifier.predict(self.X)

    def get_pred_summ(self) -> None:
        self.pred_summ = self.predictions.summarize(show_extra_info=True, verbose=True)

        ## TODO: add method to get actual class labels to use as input for `pred_summ`. Required for independent test samples
        #self.pred_summ = self.predictions.summarize(actual_classes=self.y, show_top_features=True, verbose=True)

    def get_vis_data(self):
        self.vis_data = self.predictions.get_vis_data()

    def add_filename_affixes(self, filename: str):
        if self.sample_id is not None:
            filename = f"{self.sample_id}.{filename}"

        if self.compress_tsv_files and filename.endswith(".tsv"):
            filename = f"{filename}.gz"

        return filename

    @property
    def vis_data_path(self):
        filename = self.add_filename_affixes("cuppa.vis_data.tsv")
        return os.path.join(self.output_dir, filename)

    @property
    def plot_path(self):
        filename = self.add_filename_affixes("cuppa.vis.png")
        return os.path.join(self.output_dir, filename)

    @property
    def pred_summ_path(self):
        filename = self.add_filename_affixes("cuppa.pred_summ.tsv")
        return os.path.join(self.output_dir, filename)

    def run(self) -> None:
        self.get_X()
        self.get_predictions()
        self.get_pred_summ()
        self.get_vis_data()

        self.pred_summ.to_tsv(self.pred_summ_path, verbose=True)
        self.vis_data.to_tsv(self.vis_data_path, verbose=True)

        if self.X.shape[0] > 1:
            self.logger.warning("Cannot plot visualization when predicting on multiple samples")
            return None

        plotter = CuppaVisPlotter(vis_data=self.vis_data, plot_path=self.plot_path, vis_data_path=self.vis_data_path)
        plotter.plot()