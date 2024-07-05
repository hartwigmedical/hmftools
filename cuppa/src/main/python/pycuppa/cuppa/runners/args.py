from __future__ import annotations

import argparse

from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH


class DEFAULT_RUNNER_ARGS:

    output_prefix: str = None
    compress_tsv_files: bool = False

    min_samples_with_rna: int = 5
    excl_classes: str | list[str] = ["_Other","_Unknown"]

    fusion_overrides_path: str = DEFAULT_FUSION_OVERRIDES_PATH
    clf_group: str = "all"

    skip_cv: bool = False
    cv_folds: int = 10
    cache_training: bool = True
    n_jobs: int = 1

    log_to_file: bool = False
    log_path: str | None = None
    log_format: str = "%(asctime)s [%(threadName)s] %(levelname)s %(name)s | %(message)s"

def comma_sep_str_to_list(string) -> list[str]:
    return [s.strip() for s in string.split(",")]


class RunnerArgs:

    """This class contains the kwargs to pass to argparse.ArgumentParser.add_argument()"""

    ## Required ================================
    features_path = dict(
        required=True,
        help="[Old features format] Directory with the input features csv files."
             "[New features format] Path to the input features tsv file"
    )

    output_dir = dict(
        required=True,
        help="Output directory"
    )

    metadata_path = dict(
        required=True,
        help="Path to the metadata file with cancer type labels per sample"
    )

    classifier_path = dict(
        required=True,
        help="Path to the pickle file of the classifier (e.g. cuppa_classifier.pickle.gz)"
    )

    sample_id = dict(
        help="Sample ID. Prepended to output filenames and used to look up sample in cross-validation predictions"
    )

    ## Options ================================
    clf_group = dict(
        type=str,
        default=DEFAULT_RUNNER_ARGS.clf_group,
        help=f"Classifier groups to filter probabilities for. Can be 'all' or 'dna'. Default: {DEFAULT_RUNNER_ARGS.clf_group}"
    )

    compress_tsv_files = dict(
        dest="compress_tsv_files",
        action="store_true",
        help="Compress tsv files with gzip? (will add .gz to the file extension)"
    )

    min_samples_with_rna = dict(
        type=int,
        default=DEFAULT_RUNNER_ARGS.min_samples_with_rna,
        help="Minimum number of samples with RNA in each cancer subtype. "
             "If the cancer subtype has fewer samples with RNA than this value, the cancer subtype will be excluded "
             "from training. Default: " + str(DEFAULT_RUNNER_ARGS.min_samples_with_rna)
    )

    excl_classes = dict(
        type=comma_sep_str_to_list,
        default=DEFAULT_RUNNER_ARGS.excl_classes,
        help="Comma separated list of cancer subtypes to exclude from training. "
             "E.g. 'Breast' or 'Breast,Lung'. Default: "
             ",".join(DEFAULT_RUNNER_ARGS.excl_classes)
    )

    fusion_overrides_path = dict(
        type=str,
        default=DEFAULT_RUNNER_ARGS.fusion_overrides_path,
        help="Path to the fusion overrides tsv file"
    )

    ## Cross-validation / training ================================
    cv_predictions_path = dict(
        help="Path to a CuppaPrediction tsv file containing the cross-validation predictions."
             "If provided, sample ids found in this file will have their predictions returned from this file"
             "instead of being computed"
    )

    cv_folds = dict(
        type=int,
        default=DEFAULT_RUNNER_ARGS.cv_folds,
        help="Number of cross-validation folds. Default: " + str(DEFAULT_RUNNER_ARGS.cv_folds)
    )

    skip_cv = dict(
        action="store_true",
        help="Skip cross-validation?"
    )

    no_cache_training = dict(
        dest="cache_training",
        action="store_false",
        help="Don't cache models during cross-validation and training of final model?"
    )

    n_jobs = dict(
        type=int,
        default=DEFAULT_RUNNER_ARGS.n_jobs,
        help="Number of threads to use during cross-validation (parallelized cross-validation folds) "
             "and final model training (ColumnTransformer -> parallelized sub-/meta-classifier training). "
             f"If -1, all threads are used. Default: {DEFAULT_RUNNER_ARGS.n_jobs}"
    )

    ## Logging ================================
    log_to_file = dict(
        dest="log_to_file",
        action="store_true",
        help="Output logs to a file?"
    )

    log_path = dict(
        default=DEFAULT_RUNNER_ARGS.log_path,
        help="Path to log file. Default: $(output_dir)/run.log"
    )

    log_format = dict(
        default=DEFAULT_RUNNER_ARGS.log_format,
        help="Log format"
    )


class RunnerArgParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser()

    def add_argument(self, name: str) -> None:
        kwargs = getattr(RunnerArgs, name)
        self.parser.add_argument(f"--{name}", **kwargs)

    def get_kwargs(self) -> dict:
        return self.parser.parse_args().__dict__

    def get_kwargs_predict(self) -> dict:
        self.parser.description = "PredictionRunner"

        arg_names = [
            "classifier_path",
            "features_path",
            "output_dir",
            "sample_id",
            "clf_group",
            "cv_predictions_path",
            "compress_tsv_files",
            "log_to_file",
            "log_path",
            "log_format"
        ]

        for name in arg_names:
            self.add_argument(name)

        return self.get_kwargs()

    def get_kwargs_train(self) -> dict:
        self.parser.description = "TrainingRunner"

        arg_names = [
            "features_path",
            "output_dir",
            "metadata_path",

            "min_samples_with_rna",
            "excl_classes",
            "fusion_overrides_path",

            "skip_cv",
            "cv_folds",
            "no_cache_training",
            "n_jobs",

            "log_to_file",
            "log_path",
            "log_format"
        ]

        for name in arg_names:
            self.add_argument(name)

        return self.get_kwargs()

