from __future__ import annotations

import os
from functools import cached_property

from sklearn.feature_selection import chi2

import pandas as pd
import numpy as np

from cuppa.constants import MOCK_DATA_DIR
from cuppa.logger import LoggerMixin
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, CuppaFeatures
from cuppa.sample_data.sample_metadata import SampleMetadata
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredSummary, PerformanceStats
from cuppa.compose.pipeline import PipelineCrossValidator


if __name__ == '__main__':

    from cuppa.tests.mock_data_builder import (
        MockTrainingDataBuilder,
        MockCvOutputBuilder,
        MockClassifierBuilder,
    )

    ## Input paths --------------------------------
    input_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/"
    metadata_path = input_dir + "cup_ref_sample_data.csv"
    features_dir = input_dir

    ## Output paths --------------------------------
    mock_features_path = os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz")
    mock_metadata_path = os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv")
    mock_classifier_path = os.path.join(MOCK_DATA_DIR, "cuppa_classifier.mock.pickle.gz")

    ## Training data --------------------------------
    n_samples_per_class = pd.DataFrame.from_records([
        {"class": "Breast: Other", "class_name": "Breast", "n_dna_rna": 26, "n_dna_only": 2},
        {"class": "Prostate", "class_name": "Prostate", "n_dna_rna": 22, "n_dna_only": 2},
        {"class": "Skin: Melanoma", "class_name": "Melanoma", "n_dna_rna": 18, "n_dna_only": 2},
        {"class": "Lung: Non-small cell: LUAD", "class_name": "Lung", "n_dna_rna": 12, "n_dna_only": 2},
        {"class": "Myeloid: Acute myeloid leukemia", "class_name": "AML", "n_dna_rna": 10, "n_dna_only": None}
    ])

    training_data_builder = MockTrainingDataBuilder(
        metadata_path = input_dir + "cup_ref_sample_data.csv",
        features_dir = input_dir,
        n_samples_per_class=n_samples_per_class
    )

    training_data_builder.build_and_export(os.path.join(MOCK_DATA_DIR, "training_data"))

    ## CV --------------------------------
    cv_data_builder = MockCvOutputBuilder.from_paths(
        features_path=mock_features_path,
        metadata_path=mock_metadata_path
    )
    cv_data_builder.build_and_export(out_dir=os.path.join(MOCK_DATA_DIR, "cv_output"))

    ## Classifier --------------------------------
    classifier_builder = MockClassifierBuilder.from_paths(
        features_path=mock_features_path,
        metadata_path=mock_metadata_path
    )
    classifier_builder.build_and_export(mock_classifier_path)

    # ## Intermediate probs --------------------------------
    # probs_builder = MockIntermediateProbsBuilder(
    #     cuppa_classifier_path=mock_classifier_path,
    #     features_path=mock_features_path
    # )
    # probs_builder.build_and_export(os.path.join(MOCK_DATA_DIR, "probs_fit_transform"))


class MockTrainingDataBuilder(LoggerMixin):

    def __init__(
        self,
        metadata_path: str,
        features_dir: str,
        n_samples_per_class: pd.DataFrame,
        seed: int = 0
    ):
        self.metadata_path = metadata_path
        self.features_dir = features_dir

        self._check_n_samples_per_class(n_samples_per_class)
        self.n_samples_per_class = n_samples_per_class

        self.seed = seed

    def _check_n_samples_per_class(self, df):
        required_columns = ("class", "n_dna_rna", "n_dna_only", "class_name")
        if not df.columns.isin(required_columns).all():
            self.logger.error("`n_samples_per_class` must have the columns: " + ", ".join(required_columns))
            raise LookupError

    def load_metadata(self) -> SampleMetadata:
        self.logger.info("Loading metadata from: " + self.metadata_path)
        metadata = SampleMetadata.from_csv(self.metadata_path)
        metadata["data_type"] = (metadata["RnaReadLength"] != 0).map({True: "dna_rna", False: "dna_only"})
        return metadata

    def select_samples_by_data_type(self, metadata: SampleMetadata,  data_type: str) -> SampleMetadata:
        #data_type = "dna_rna"

        self.logger.info(f"Selecting '{data_type}' samples")

        sample_counts = metadata[["CancerSubtype", "data_type"]] \
            .value_counts() \
            .unstack() \
            .fillna(0) \
            .astype(int)

        valid_data_types = ("dna_rna", "dna_only")
        if data_type not in valid_data_types:
            self.logger.error("`data_type` must be one of: " + ", ".join(valid_data_types))
            raise ValueError

        rng = np.random.default_rng(self.seed)
        selected_samples = {}

        for index, row in self.n_samples_per_class.iterrows():
            # break
            class_i = row["class"]
            sample_ids_i = metadata.index[
                (metadata["CancerSubtype"] == class_i)
                & (metadata["data_type"] == data_type)
            ]

            actual_n_samples_i = sample_counts.loc[class_i, data_type]
            target_n_samples_i = row[f"n_{data_type}"]

            if np.isnan(target_n_samples_i) or target_n_samples_i >= actual_n_samples_i:
                selected_samples_i = pd.Series(sample_ids_i)
            else:
                target_n_samples_i = int(target_n_samples_i)
                selected_indexes = rng.integers(low=0, high=actual_n_samples_i, size=target_n_samples_i)
                selected_samples_i = pd.Series(sample_ids_i[selected_indexes])

            self.logger.info(f"Selected {len(selected_samples_i)} samples for '{class_i}' for '{data_type}'")
            selected_samples[class_i] = selected_samples_i

        selected_samples = pd.concat(selected_samples).values.tolist()

        return metadata.loc[selected_samples]

    def load_features(self) -> CuppaFeatures:
        self.logger.info("Loading features from: " + self.features_dir)

        features_paths = CuppaFeaturesPaths.from_dir(self.features_dir, basenames_mode="old")

        features_loader = FeatureLoaderOld(
            paths = features_paths,
            genome_version = 37,
            excl_chroms = ["ChrY", "Y"],
            verbose = True
        )

        return features_loader.load_features()

    def filter_na_zero_columns(self, X: CuppaFeatures) -> CuppaFeatures:
        is_na_or_zero = X.isna() | (X==0)
        return X.loc[:,~is_na_or_zero.all()]


    def filter_irrelevant_features(
        self,
        X: CuppaFeatures,
        y: pd.Series,
        pvalue_thres: float,
        feature_type: str
    ) -> CuppaFeatures:
        # X = features.get_feat_type_cols("alt_sj")
        # y = metadata.index
        # top_n = 20_000

        chi2_stat, pvalue = chi2(X.fillna(0), y)

        test_results = pd.DataFrame(
            dict(chi2_stat=chi2_stat, pvalue=pvalue),
            index=X.columns
        )

        test_results["chi2_stat"] = test_results["chi2_stat"].fillna(0)
        test_results["pvalue"] = test_results["pvalue"].fillna(1)

        sel_features = test_results.index[test_results["pvalue"] < pvalue_thres]
        sel_features = X.columns[X.columns.isin(sel_features)]

        self.logger.info("Selected %i/%i features for %s" % (
            len(sel_features), len(test_results),
            feature_type
        ))

        return X[sel_features]

    def format_metadata(self, metadata: SampleMetadata) -> SampleMetadata:

        metadata = metadata.copy()[["CancerType", "CancerSubtype", "RnaReadLength"]]

        self.logger.info("Renaming classes")
        class_mappings = dict(zip(
            self.n_samples_per_class["class"],
            self.n_samples_per_class["class_name"]
        ))
        metadata["CancerSubtype"] = metadata["CancerSubtype"].map(class_mappings)

        self.logger.info("Renaming samples")
        metadata.index = \
            pd.Series(range(0, len(metadata))).astype(str) \
            + "_" \
            + metadata["CancerSubtype"].values.astype(str)

        return metadata

    @cached_property
    def _build_cached(self) -> tuple[SampleMetadata, CuppaFeatures]:

        ## --------------------------------
        metadata_full = self.load_metadata()

        metadata = {}
        metadata["dna_rna"] = self.select_samples_by_data_type(metadata_full, "dna_rna")
        metadata["dna_only"] = self.select_samples_by_data_type(metadata_full, "dna_only")

        metadata = pd.concat(metadata.values())

        ## --------------------------------
        features_full = self.load_features()

        features = features_full.loc[metadata.index]

        self.logger.info("Removing feature columns with all NAs or zeros")
        features = self.filter_na_zero_columns(features)

        is_rna_feature = features.columns.str.match("^(gene_exp|alt_sj)")
        features_dna = features.loc[:,~is_rna_feature]
        features_rna = {}

        features_rna["gene_exp"] = self.filter_irrelevant_features(
            X=features.loc[:, features.columns.str.startswith("gene_exp")],
            y=metadata.index,
            pvalue_thres=0.05,
            feature_type="gene_exp"
        )

        features_rna["alt_sj"] = self.filter_irrelevant_features(
            X=features.loc[:,features.columns.str.startswith("alt_sj")],
            y=metadata.index,
            pvalue_thres=0.05,
            feature_type="alt_sj"
        )

        features_rna = pd.concat(features_rna.values(), axis=1)

        features = pd.concat([features_dna, features_rna], axis=1)

        ## --------------------------------
        metadata = self.format_metadata(metadata)

        features.index = metadata.index

        return metadata, features

    def build_and_export(self, out_dir: str) -> None:

        metadata = self._build_cached[0]
        features = self._build_cached[1]

        self.logger.info("Exporting output to: " + out_dir)
        os.makedirs(out_dir, exist_ok=True)
        metadata.to_csv(os.path.join(out_dir, "metadata.tsv"), sep='\t')
        features.to_csv(os.path.join(out_dir, "features.tsv.gz"), sep='\t')


class MockCvOutputBuilder(LoggerMixin):
    def __init__(self, features: CuppaFeatures, metadata: SampleMetadata):
        self.X = features
        self.y = metadata["CancerSubtype"]
        self.y_split = metadata["CancerSubtype"] + "__" + metadata["RnaReadLength"].astype(str)

    @classmethod
    def from_paths(cls, features_path: str, metadata_path: str) -> "MockCvOutputBuilder":
        features = pd.read_csv(features_path, sep='\t', index_col=0)
        metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
        return cls(features=features, metadata=metadata)

    def cross_validate(self) -> PipelineCrossValidator:

        cross_validator = PipelineCrossValidator(
            pipeline=CuppaClassifier(fusion_overrides_path=None),
            X=self.X,
            y=self.y,
            y_split=self.y_split,
        )
        cross_validator.fit(cache_training=False)
        return cross_validator

    @cached_property
    def _build_cached(self) -> tuple[CuppaPrediction, CuppaPredSummary, PerformanceStats]:
        cross_validator = self.cross_validate()
        #probs = cross_validator.apply_on_test_sets("predict_proba", n_jobs=1)
        predictions = cross_validator.apply_on_test_sets("predict", probs_only=False, n_jobs=1)
        pred_summ = predictions.summarize(actual_classes=self.y)
        performance = pred_summ.performance()

        return (predictions, pred_summ, performance)

    def build_and_export(self, out_dir: str):
        predictions, pred_summ, performance = self._build_cached

        self.logger.info("Exporting output to: " + out_dir)
        os.makedirs(out_dir, exist_ok=True)

        predictions.to_tsv(os.path.join(out_dir, "predictions.tsv.gz"))
        pred_summ.to_tsv(os.path.join(out_dir, "pred_summ.tsv"))
        performance.to_tsv(os.path.join(out_dir, "performance.tsv"))


class MockClassifierBuilder:
    def __init__(self, features: CuppaFeatures, metadata: SampleMetadata):
        self.X = features
        self.y = metadata["CancerSubtype"]

    @classmethod
    def from_paths(cls, features_path: str, metadata_path: str) -> "MockClassifierBuilder":
        features = pd.read_csv(features_path, sep='\t', index_col=0)
        metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
        return cls(features=features, metadata=metadata)

    def build_and_export(self, out_path: str):
        classifier = CuppaClassifier(fusion_overrides_path=None)
        classifier.fit(X=self.X, y=self.y)
        classifier.to_file(out_path)

#
# class MockIntermediateProbsBuilder(LoggerMixin):
#     def __init__(self, cuppa_classifier_path: str, features_path: str):
#         self.cuppa_classifier = CuppaClassifier.from_file(cuppa_classifier_path)
#         self.features = pd.read_csv(features_path, sep='\t', index_col=0)
#
#     @cached_property
#     def _build_cached(self) -> dict[str, pd.DataFrame]:
#
#         probs = self.cuppa_classifier.predict_proba(self.features)
#
#         clf_names = probs.index.get_level_values("clf_name").unique()
#         probs_per_clf = {}
#         for clf_name in clf_names:
#             probs_i = probs[probs.index.get_level_values("clf_name") == clf_name].copy()
#             probs_i.index = probs_i.index.get_level_values("sample_id")
#             probs_i.columns.name = None
#             probs_per_clf[clf_name] = probs_i
#
#         return probs_per_clf
#
#     def build_and_export(self, out_dir: str):
#         #out_dir = os.path.join(MOCK_DATA_DIR, "probs_fit_transform/")
#
#         self.logger.info("Exporting output to: " + out_dir)
#         os.makedirs(out_dir, exist_ok=True)
#
#         probs_per_clf = self._build_cached
#         for clf_name, probs_i in probs_per_clf.items():
#             probs_i.to_csv(os.path.join(out_dir, clf_name + ".tsv"), sep='\t')