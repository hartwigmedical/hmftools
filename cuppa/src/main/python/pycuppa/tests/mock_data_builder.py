from __future__ import annotations

import os
from functools import cached_property

import pandas as pd
import numpy as np
from sklearn.feature_selection import chi2
from sklearn.model_selection import StratifiedKFold

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
        MockOldCuppaFeaturesBuilder
    )

    ## Input paths --------------------------------
    input_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/"
    metadata_path = input_dir + "cup_ref_sample_data.csv"
    features_dir = input_dir

    ## Output paths --------------------------------
    mock_features_path = os.path.join(MOCK_DATA_DIR, "training_data/features.tsv.gz")
    mock_metadata_path = os.path.join(MOCK_DATA_DIR, "training_data/metadata.tsv")

    ## Training data --------------------------------
    n_samples_per_class = pd.DataFrame.from_records([
        {"class": "Breast: Other", "class_name": "Breast", "n_dna_rna": 26, "n_dna_only": 2},
        {"class": "Prostate", "class_name": "Prostate", "n_dna_rna": 22, "n_dna_only": 2},
        {"class": "Skin: Melanoma", "class_name": "Melanoma", "n_dna_rna": 18, "n_dna_only": 2},
        {"class": "Lung: Non-small cell: LUAD", "class_name": "Lung", "n_dna_rna": 12, "n_dna_only": 2},
        {"class": "Myeloid: Acute myeloid leukemia", "class_name": "AML", "n_dna_rna": 10, "n_dna_only": None}
    ])

    training_data_builder = MockTrainingDataBuilder(
        metadata_path = metadata_path,
        features_dir = input_dir,
        n_samples_per_class=n_samples_per_class
    )

    training_data_builder.build_and_export(os.path.join(MOCK_DATA_DIR, "training_data"))

    fusion_overrides = pd.DataFrame.from_records(
        columns=("feat_prefix", "feat_basename", "target_class"),
        data=[
            ("event.fusion.", "CBFB_MYH11", "AML"),
            ("event.fusion.", "RUNX1_RUNX1T1", "AML"),
        ]
    )
    fusion_overrides.to_csv(
        os.path.join(MOCK_DATA_DIR, "training_data/fusion_overrides.tsv"),
        sep='\t', index=False
    )

    ## CV --------------------------------
    cv_data_builder = MockCvOutputBuilder.from_paths(
        features_path=mock_features_path,
        metadata_path=mock_metadata_path
    )
    cv_data_builder.build_and_export(out_dir=os.path.join(MOCK_DATA_DIR, "cv_output"))

    ## Final model --------------------------------
    classifier_builder = MockClassifierBuilder.from_paths(
        features_path=mock_features_path,
        metadata_path=mock_metadata_path
    )
    classifier_builder.build_and_export(out_dir=os.path.join(MOCK_DATA_DIR, "final_model"))

    ## --------------------------------
    mock_features_builder = MockOldCuppaFeaturesBuilder(out_dir=os.path.join(MOCK_DATA_DIR, "input_data/old_format"))
    mock_features_builder.build_and_export()


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

        features_paths = CuppaFeaturesPaths.from_dir(self.features_dir, file_format="old")

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
            cv=StratifiedKFold(n_splits=10, random_state=0, shuffle=True),
        )
        cross_validator.fit(cache_training=False)
        return cross_validator

    @cached_property
    def _build_cached(self) -> tuple[CuppaPrediction, CuppaPredSummary, PerformanceStats]:
        cross_validator = self.cross_validate()
        #probs = cross_validator.apply_on_test_sets("predict_proba", n_jobs=1)
        predictions = cross_validator.apply_on_test_sets("predict", n_jobs=1)
        pred_summ = predictions.summarize(actual_classes=self.y)
        performance = pred_summ.performance()

        return (predictions, pred_summ, performance)

    def build_and_export(self, out_dir: str) -> None:
        predictions, pred_summ, performance = self._build_cached

        self.logger.info("Exporting output to: " + out_dir)
        os.makedirs(out_dir, exist_ok=True)

        predictions.to_tsv(os.path.join(out_dir, "predictions.tsv.gz"))
        pred_summ.to_tsv(os.path.join(out_dir, "pred_summ.tsv"))
        performance.to_tsv(os.path.join(out_dir, "performance.tsv"))


class MockClassifierBuilder(LoggerMixin):
    def __init__(self, features: CuppaFeatures, metadata: SampleMetadata):
        self.X = features
        self.y = metadata["CancerSubtype"]

    @classmethod
    def from_paths(cls, features_path: str, metadata_path: str) -> "MockClassifierBuilder":
        features = pd.read_csv(features_path, sep='\t', index_col=0)
        metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
        return cls(features=features, metadata=metadata)

    def build_and_export(self, out_dir: str) -> None:
        classifier = CuppaClassifier(fusion_overrides_path=None)
        classifier.fit(X=self.X, y=self.y)

        predictions = classifier.predict(self.X)
        predictions.get_data_types("prob").get_clf_names("gen_pos")

        self.logger.info("Exporting output to: " + out_dir)
        os.makedirs(out_dir, exist_ok=True)

        predictions.to_tsv(os.path.join(out_dir, "predictions.tsv.gz"), signif_digits=None)
        classifier.to_file(os.path.join(out_dir, "cuppa_classifier.pickle.gz"))


class MockOldCuppaFeaturesBuilder(LoggerMixin):

    def __init__(self, out_dir: str):
        self.out_dir = out_dir

    def get_out_path(self, feature_type: str):
        BASENAMES = dict(
            gen_pos="cup_ref_sample_pos_freq_counts.csv.gz",
            snv96="cup_ref_snv_counts.csv.gz",
            driver_fusion_virus="cup_ref_cohort_feature_data.csv",
            sv="cup_ref_cohort_sv_data.csv",
            trait="cup_ref_cohort_traits_data.csv",
            sig="cup_ref_cohort_signature_data.csv",

            gene_exp="cup_ref_gene_exp_sample.csv",
            alt_sj="cup_ref_alt_sj_sample.csv",
        )
        return os.path.join(self.out_dir, BASENAMES[feature_type])


    @staticmethod
    def make_dummy_names(n: int, prefix: str = "sample_"):
        return prefix + pd.Series(range(1, n+1)).astype(str)

    def make_feature_x_sample_matrix(
        self,
        n_samples: int = 2,
        n_features: int = 3,
        dtype: object = int,
        feature_names: list[str] | pd.Index | None = None
    ) -> pd.DataFrame:
        #n_samples=3
        #n_features=2

        matrix = pd.DataFrame(
            np.zeros((n_samples, n_features)),
            index = "sample_" + pd.Series(range(1, n_samples+1)).astype(str)
        )

        if feature_names is not None:
            matrix.columns = feature_names

        matrix = matrix.astype(dtype)

        return matrix

    def make_gen_pos_matrix(self, n_samples: int = 2) -> pd.DataFrame:
        return self.make_feature_x_sample_matrix(n_samples=n_samples, n_features=6220).transpose()

    def make_snv96_matrix(self, n_samples: int = 2) -> pd.DataFrame:
        return self.make_feature_x_sample_matrix(n_samples=n_samples, n_features=96).transpose()

    def make_sv_matrix(self, n_samples: int = 2) -> pd.DataFrame:
        matrix = self.make_feature_x_sample_matrix(
            n_samples=n_samples, n_features=6,
            feature_names=[
                "LINE",
                "SIMPLE_DEL_20KB_1MB","SIMPLE_DUP_32B_200B","SIMPLE_DUP_100KB_5MB",
                "MAX_COMPLEX_SIZE","TELOMERIC_SGL"
            ]
        )
        matrix = matrix.reset_index(names="SampleId")
        return matrix

    def make_trait_matrix(self, n_samples: int = 2) -> pd.DataFrame:
        matrix = self.make_feature_x_sample_matrix(
            n_samples=n_samples, n_features=6, dtype=float,
            feature_names=["Gender", "WholeGenomeDuplication", "Purity", "Ploidy", "MsIndelsPerMb", "ChordHrd"]
        )
        matrix = matrix.reset_index(names="SampleId")

        matrix["Gender"] = "FEMALE"
        matrix["WholeGenomeDuplication"] = "false"

        return matrix

    def make_gene_exp_matrix(self, n_samples: int = 2) -> pd.DataFrame:

        feature_names = pd.MultiIndex.from_tuples([
            ("ENSG00000000001", "GENE1"),
            ("ENSG00000000002", "GENE2"),
            ("ENSG00000000003", "GENE3"),
        ], names=["GeneId","GeneName"])

        matrix = self.make_feature_x_sample_matrix(
            n_samples=n_samples, n_features=len(feature_names), dtype=float,
            feature_names=feature_names
        )

        matrix = matrix.transpose().reset_index()

        return matrix

    def make_alt_sj_matrix(self, n_samples: int = 2) -> pd.DataFrame:

        feature_names = pd.MultiIndex.from_tuples([
            ("ENSG00000000001", "1", 1000, 2000),
            ("ENSG00000000002", "1", 3000, 4000),
            ("ENSG00000000003", "1", 5000, 6000),
        ], names=["GeneId", "Chromosome","PosStart","PosEnd"])

        matrix = self.make_feature_x_sample_matrix(
            n_samples=n_samples, n_features=len(feature_names),
            feature_names=feature_names
        )

        matrix = matrix.transpose().reset_index()

        return matrix

    def make_event_matrix_one_sample(
        self,
        event_types = ["DRIVER", "FUSION", "VIRUS", "INDEL"],
        sample_num: int = 1
    ) -> pd.DataFrame:

        template = pd.DataFrame.from_dict(dict(
            SampleId="sample_" + str(sample_num),
            Name=["driver_1", "fusion_1", "virus_1", "indel_1"],
            Type=["DRIVER", "FUSION", "VIRUS", "INDEL"],
            Likelihood=[0.0, 0.0, 0.0, 0.0],
            ExtraInfo=["TYPE=DEL", "", "", ""]
        ))

        return template[template["Type"].isin(event_types)]

    def make_signature_matrix(self, n_samples: int = 2) -> pd.DataFrame:
        template = pd.DataFrame.from_dict(dict(
            SampleId="sample_",
            SigName=["Sig2", "Sig13", "Sig4", "Sig6", "Sig7", "Sig17"],
            Allocation=0.0
        ))

        dfs = []
        for i in range(1, n_samples + 1):
            df = template.copy()
            df["SampleId"] = df["SampleId"] + str(i)
            dfs.append(df)

        return pd.concat(dfs).reset_index(drop=True)

    def build_and_export(self):

        n_dna_samples = 3
        self.make_gen_pos_matrix(n_samples=n_dna_samples).to_csv(self.get_out_path("gen_pos"), index=False)
        self.make_snv96_matrix(n_samples=n_dna_samples).to_csv(self.get_out_path("snv96"), index=False)
        self.make_sv_matrix(n_samples=n_dna_samples).to_csv(self.get_out_path("sv"), index=False)
        self.make_signature_matrix(n_samples=n_dna_samples).to_csv(self.get_out_path("sig"), index=False)
        self.make_trait_matrix(n_samples=n_dna_samples).to_csv(self.get_out_path("trait"), index=False)

        self.make_event_matrix_one_sample().to_csv(self.get_out_path("driver_fusion_virus"), index=False)

        n_rna_samples = 2
        self.make_gene_exp_matrix(n_samples=n_rna_samples).to_csv(self.get_out_path("gene_exp"), index=False)
        self.make_alt_sj_matrix(n_samples=n_rna_samples).to_csv(self.get_out_path("alt_sj"), index=False)




