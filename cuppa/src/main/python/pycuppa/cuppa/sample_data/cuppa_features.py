from __future__ import annotations

import importlib
import os.path
from functools import cached_property
from typing import Iterable, Literal, Self

import numpy as np
import pandas as pd

from cuppa.components.preprocessing import NaRowFilter
from cuppa.logger import LoggerMixin


class CuppaFeaturesPaths(pd.Series, LoggerMixin):

    def __init__(self, paths: dict[str, str] | pd.Series):
        super().__init__(paths)

    BASENAMES_OLD = dict(
        #metadata="cup_ref_sample_data.csv",
        gen_pos="cup_ref_sample_pos_freq_counts.csv",
        snv96="cup_ref_snv_counts.csv",
        driver_fusion_virus="cup_ref_cohort_feature_data.csv",
        sv="cup_ref_cohort_sv_data.csv",
        trait="cup_ref_cohort_traits_data.csv",
        sig="cup_ref_cohort_signature_data.csv",

        gene_exp="cup_ref_gene_exp_sample.csv",
        alt_sj="cup_ref_alt_sj_sample.csv",
    )

    BASENAMES_NEW = dict(
        gen_pos="gen_pos.tsv",
        snv96="snv96.tsv",
        event="event.tsv",
        sig="sig.tsv",

        gene_exp="gene_exp.tsv",
        alt_sj="alt_sj.tsv",
    )

    OPTIONAL_BASENAME_KEYS = ("gene_exp", "alt_sj")

    @classmethod
    def from_dir(
        cls,
        directory: str,
        basenames_mode: Literal["old", "new"] = "new"
    ):

        if basenames_mode == "new":
            basenames_expected = cls.BASENAMES_NEW
        elif basenames_mode == "old":
            basenames_expected = cls.BASENAMES_OLD
        else:
            cls.get_class_logger(cls).error("`basenames_mode` must be 'old' or 'new'")
            raise ValueError

        basenames = {}
        basenames_missing = {}
        paths = {}

        for key, basename in basenames_expected.items():

            basename_gz = basename + ".gz"
            path = os.path.join(directory, basename_gz)
            if os.path.exists(path):
                basenames[key] = basename_gz
                paths[key] = path
                continue

            path = os.path.join(directory, basename)
            if os.path.exists(path):
                basenames[key] = basename
                paths[key] = path
                continue

            basenames_missing[key] = basename

        if len(basenames_missing)==0:
            return CuppaFeaturesPaths(paths)


        required_basenames_missing = {}
        logger = cls.get_class_logger(cls)
        for key, basename in basenames_missing.items():

            if key in cls.OPTIONAL_BASENAME_KEYS:
                logger.warning("Missing optional input file for %s: %s" % (key, basename))
                continue

            logger.error("Missing required input file for %s: %s" % (key, basename))
            required_basenames_missing[key] = basename

        if len(required_basenames_missing) > 0:
            logger.info("Could not load paths to features. Maybe set basenames_mode='old'?")
            raise Exception

        return CuppaFeaturesPaths(paths)


class CuppaFeatures(pd.DataFrame, LoggerMixin):
    def __init__(self, df: pd.DataFrame, *args, **kwargs):
        super().__init__(df, *args, **kwargs)

    @property
    def _constructor(self):
        return CuppaFeatures

    @cached_property
    def col_feat_types(self) -> pd.Index:
        ## Match up to but excluding first dot
        return self.columns.str.extract("(^\w+[^.])", expand=False)

    @cached_property
    def feat_types(self) -> pd.Series:
        return pd.Series(self.col_feat_types.unique())

    def get_feat_type_cols(self, sel_feat_types: str | Iterable[str]):
        sel_feat_types = pd.Series(sel_feat_types)

        invalid_feat_types = sel_feat_types[~sel_feat_types.isin(self.feat_types)]
        if len(invalid_feat_types)>0:
            self.logger.error("Invalid feature types: " + ", ".join(invalid_feat_types))
            raise LookupError

        sel_cols = self.col_feat_types.isin(sel_feat_types)
        return self.loc[:, sel_cols]

    FEAT_TYPES_WIDE = ["gen_pos", "gene_exp", "alt_sj"]

    def to_tsv_files(
        self,
        out_dir: str,
        drop_na_rows: bool = True,
        float_format: str = "%.8g",
        chunksize: int = 1000,
        verbose: bool = True,
        *args, **kwargs
    ) -> None:
        # out_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/output/cuppa_features/"

        for feat_type in self.feat_types:
            # feat_type="gene_exp"

            path = os.path.join(out_dir, feat_type + ".tsv.gz")
            df = self.get_feat_type_cols(feat_type)

            if drop_na_rows:
                na_rows = NaRowFilter.detect_na_rows(df, use_first_col=True)
                df = df.loc[~na_rows]

            if feat_type in self.FEAT_TYPES_WIDE:
                df = df.transpose()

            if verbose:
                self.logger.info("Writing %s tsv to: %s" % (feat_type, path))

            df.to_csv(
                path,
                sep="\t",
                float_format=float_format,
                chunksize=chunksize,
                *args, **kwargs
            )

    # TODO: the below method is subject to change
    @classmethod
    def from_tsv_files(cls, paths: CuppaFeaturesPaths, verbose: bool = True) -> Self:

        features = {}
        for feat_type, path in paths.items():
            if verbose:
                cls.get_class_logger(cls).info("Reading %s tsv from: %s" % (feat_type, path))

            df = pd.read_table(path, index_col=0, engine="c")
            if feat_type in cls.FEAT_TYPES_WIDE:
                df = df.transpose()

            features[feat_type] = df

        features = pd.concat(features.values(), axis=1)
        return CuppaFeatures(features)


class FeatureLoaderNew(LoggerMixin):

    def __init__(self, path: str, verbose: bool = False):
        self.path = path
        self.verbose = verbose

    ## Loading ================================
    FEAT_INFO_COLS = dict(
        Source = "source",
        Category = "category",
        Key = "key"
    )

    def load_feat_info(self, path: str) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Loading features info")

        df = pd.read_table(
            path,
            index_col=False,
            usecols=list(self.FEAT_INFO_COLS.keys())
        )

        df.columns = self.FEAT_INFO_COLS.values()

        return df

    ## Helper methods to parse feature names ================================
    def _rename_categories(self, feat_info: pd.DataFrame) -> None:

        if self.verbose:
            self.logger.debug("Renaming `Category` values")

        mappings = dict(
            GEN_POS="gen_pos",
            SNV96="snv96",
            SV_COUNT="event.sv",
            FUSION="event.fusion",
            SAMPLE_TRAIT="event.trait",
            SIGNATURE="sig",

            EXPRESSION="gene_exp",
            ALT_SJ="alt_sj",

            ## TODO: below mappings are subject to change
            DRIVER="event.driver",
            VIRUS="virus",
        )

        feat_info["category_renamed"] = feat_info["category"].replace(mappings)


    def _rename_keys(self, feat_info: pd.DataFrame) -> None:

        if self.verbose:
            self.logger.debug("Renaming `Key` values")

        mappings_sigs = dict(
            SIG_1 = "Age (SBS1)",
            SIG_2_13_AID_APOBEC = "AID/APOBEC (SBS2/13)",
            SIG_4_SMOKING = "Smoking (SBS4)",
            SIG_6_MMR = "MMRD (SBS6)",
            SIG_7_UV = "UV (SBS7)",
            SIG_10_POLE = "POLE (SBS10)",
            SIG_11 = "Temozolomide (SBS11)",
            SIG_17 = "ROS/5FU (SBS17)",
        )

        mappings_traits = dict(
            IS_MALE = "is_male",
            SNV_COUNT="tmb.snv_count",
            MS_INDELS_TMB = "tmb.indels_per_mb",
            WGD = "whole_genome_duplication",
        )

        mappings = mappings_sigs | mappings_traits
        feat_info["key_renamed"] = feat_info["key"].replace(mappings)

    def _make_final_feat_names(self, feat_info: pd.DataFrame) -> None:

        if self.verbose:
            self.logger.debug("Renaming `feat_name` values")

        mappings = {
            "event.trait.tmb.snv_count": "event.tmb.snv_count",
            "event.trait.tmb.indels_per_mb": "event.tmb.indels_per_mb"
        }

        feat_names = feat_info["category_renamed"] + "." + feat_info["key_renamed"]
        feat_names.replace(mappings, inplace=True)

        feat_info["feat_name"] = feat_names

    def _get_feat_types(self, feat_info: pd.DataFrame, sep=".") -> None:
        affixes = feat_info["feat_name"].str.split(sep, regex=False, n=1)
        prefixes = affixes.map(lambda x: x[0])
        feat_info["feat_type"] = prefixes

    EXCLUDED_FEATURES = (
        "Age (SBS1)",
        "POLE (SBS10)",
        "Temozolomide (SBS11)"
    )

    def _mark_excluded_features(self, feat_info: pd.DataFrame) -> None:
        if self.verbose:
            self.logger.debug("Excluding %i features: %s" % (
                len(self.EXCLUDED_FEATURES),
                ", ".join(self.EXCLUDED_FEATURES)
            ))

        feat_info["is_excluded"] = feat_info["key_renamed"]\
            .isin(self.EXCLUDED_FEATURES)\
            .values

    def _mark_duplicate_features(self, feat_info: pd.DataFrame):

        feat_info["is_duplicated"] = feat_info["key_renamed"].duplicated()
        duplicate_features = feat_info.query("is_duplicated")["feat_name"]

        if self.verbose and len(duplicate_features)>0:
            self.logger.warning("Removing %i duplicate features: %s" % (
                len(duplicate_features),
                ", ".join(duplicate_features)
            ))

    ## Main ================================
    def parse_feature_names(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.info("Parsing feature names")

        feat_info = self.load_feat_info(self.path)

        self._rename_categories(feat_info)
        self._rename_keys(feat_info)
        self._make_final_feat_names(feat_info)
        self._get_feat_types(feat_info)
        self._mark_excluded_features(feat_info)
        self._mark_duplicate_features(feat_info)

        return feat_info.drop(["category_renamed", "key_renamed"], axis=1)

    def load_feature_values(self) -> pd.DataFrame:
        header = pd.read_table(self.path, nrows=0).columns
        sample_cols = header[~header.isin(self.FEAT_INFO_COLS.keys())]

        return pd.read_table(
            self.path,
            index_col=False,
            usecols=sample_cols,
            dtype="float32",
            engine='c'
        )

    def load(self):
        df = self.load_feature_values()

        ## Assign feature names
        feat_info = self.parse_feature_names()
        df.index = feat_info["feat_name"].values

        ## Remove some features
        df = df.loc[
            ~feat_info["is_excluded"].values &
            ~feat_info["is_duplicated"].values
        ]
        df = df.transpose()

        return CuppaFeatures(df)


class FeatureLoaderOld(LoggerMixin):

    def __init__(
        self,
        paths: CuppaFeaturesPaths,
        genome_version: int = 37,
        excl_chroms: Iterable[str] = ["ChrY", "Y"],
        verbose: bool = True
    ):
        self.paths = paths
        self.genome_version = genome_version
        self.excl_chroms = excl_chroms
        self.verbose = verbose

    @staticmethod
    def _add_prefix_to_columns(df, prefix: str, sep=".") -> None:
        df.columns = prefix + sep + df.columns

    ## SNV features ================================
    @cached_property
    def gen_pos_matrix(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.info("Loading features: gen_pos")

        if self.genome_version not in [37, 38]:
            self.logger.error("`genome_version`")
            raise ValueError

        feat_names = pd.read_csv(
            importlib.resources.files("resources") / ("feature_names/gen_pos.hg%i.csv" % self.genome_version)
        )

        matrix = pd.read_csv(self.paths["gen_pos"]).transpose()
        matrix.columns = feat_names["chrom"] + "_" + feat_names["pos"].astype(str)

        if self.excl_chroms is not None:
            excl_chroms = pd.Series(self.excl_chroms)
            matrix = matrix.loc[:, ~matrix.columns.str.startswith(tuple(excl_chroms))]

        return matrix

    @cached_property
    def snv96_matrix(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.info("Loading features: snv96")

        matrix = pd.read_csv(self.paths["snv96"]).transpose()

        feat_names = pd.read_csv(importlib.resources.files("resources") / "feature_names/snv96.csv")
        matrix.columns = feat_names["context"].values

        return matrix

    @cached_property
    def sig_matrix(self):

        if self.verbose:
            self.logger.info("Loading features: signatures")

        df = pd.read_csv(self.paths["sig"])

        matrix = df.pivot_table(index=df["SampleId"], columns="SigName", values="Allocation", fill_value=0)

        matrix["Sig2/Sig13"] = matrix["Sig2"] + matrix["Sig13"]

        feat_names = pd.Series({
            "Sig2/Sig13": "AID/APOBEC (SBS2/13)",
            "Sig4": "Smoking (SBS4)",
            "Sig6": "MMRD (SBS6)",
            "Sig7": "UV (SBS7)",
            "Sig17": "ROS/5FU (SBS17)"

        })

        matrix = matrix[feat_names.index]
        matrix.columns = feat_names.values

        return matrix

    ## Event features ================================
    @cached_property
    def sv_matrix(self) -> pd.DataFrame:
        matrix = pd.read_csv(self.paths["sv"], index_col="SampleId")
        self._add_prefix_to_columns(matrix, "sv")
        return matrix


    @cached_property
    def _df_driver_fusion_virus(self) -> pd.DataFrame:
        return pd.read_csv(self.paths["driver_fusion_virus"])

    @cached_property
    def driver_matrix(self) -> pd.DataFrame:

        df = self._df_driver_fusion_virus.copy()
        df = df[df["Type"].isin(["DRIVER", "INDEL"])]

        ## Remove 'INDEL' prefix
        df.loc[df["Type"] == "INDEL", "Name"] = \
            df.loc[df["Type"] == "INDEL", "Name"].str.replace("^.+_", "", regex=True)

        ## Mutation subtype
        df["sub_type"] = df["ExtraInfo"] \
            .str.extract("(TYPE=.+$)", expand=False) \
            .str.replace("TYPE=", "", regex=False) \
            .str.lower() \
            .fillna("mut")

        df.loc[df["sub_type"] == "del", "sub_type"] = "mut"
        df.loc[df["Type"] == "INDEL", "sub_type"] = "indel"

        ## Select max likelihood feature
        df = df.sort_values(["Name", "SampleId", "Likelihood"])
        df["sample_feature_id"] = df["SampleId"] + "_" + df["Name"]
        df = df.drop_duplicates("sample_feature_id", keep="last")

        ## Wide format
        df["feature_name"] = df["Name"] + "." + df["sub_type"]
        matrix = df.pivot_table(index='SampleId', columns='feature_name', values='Likelihood', fill_value=0)

        self._add_prefix_to_columns(matrix, "driver")

        return matrix

    @cached_property
    def fusion_matrix(self) -> pd.DataFrame:
        df = self._df_driver_fusion_virus.query("Type=='FUSION'")
        matrix = df.pivot_table(index='SampleId', columns='Name', values='Likelihood', fill_value=0)
        self._add_prefix_to_columns(matrix, "fusion")
        return matrix

    @cached_property
    def virus_matrix(self) -> pd.DataFrame:
        df = self._df_driver_fusion_virus.query("Type=='VIRUS'")
        matrix = df.pivot_table(index='SampleId', columns='Name', values='Likelihood', fill_value=0)
        self._add_prefix_to_columns(matrix, "virus")
        return matrix


    @cached_property
    def _df_trait(self):
        return pd.read_csv(self.paths["trait"], index_col="SampleId")

    @cached_property
    def tmb_matrix(self):
        matrix = pd.DataFrame({
            "snv_count": self.snv96_matrix.sum(axis=1),
            "indels_per_mb": self._df_trait["MsIndelsPerMb"]
        })
        self._add_prefix_to_columns(matrix, "tmb")
        return matrix

    @cached_property
    def trait_matrix(self):
        matrix = pd.DataFrame({
            "is_male": (self._df_trait["Gender"] == "MALE").astype(int),
            "whole_genome_duplication": self._df_trait["WholeGenomeDuplication"].astype(int)
        })
        self._add_prefix_to_columns(matrix, "trait")
        return matrix


    @cached_property
    def event_matrix(self):

        if self.verbose:
            self.logger.info("Loading features: event")

        l = [
            self.driver_matrix,
            self.fusion_matrix,
            self.virus_matrix,
            self.sv_matrix,
            self.tmb_matrix,
            self.trait_matrix,
        ]

        return pd.concat(l, axis=1)

    @cached_property
    def gene_exp_matrix(self):

        if self.verbose:
            self.logger.info("Loading features: gene_exp")

        path = self.paths["gene_exp"]

        ## Get sample columns
        index_cols = ["GeneId", "GeneName"]
        cols = pd.read_csv(path, nrows=0).columns
        sample_cols = cols[~cols.isin(index_cols)]

        ## Load feature values
        matrix = pd.read_csv(
            path,
            low_memory=False,
            usecols=list(sample_cols) + ["GeneName"],
            index_col="GeneName",
            engine="c"
        )
        matrix = matrix.transpose()
        matrix.fillna(0, inplace=True)

        self._add_prefix_to_columns(matrix, "gene_exp")

        return matrix

    @cached_property
    def alt_sj_matrix(self):

        if self.verbose:
            self.logger.info("Loading features: alt_sj")

        path = self.paths["alt_sj"]

        ## Get sample columns
        index_cols = ["GeneId", "Chromosome", "PosStart", "PosEnd"]
        cols = pd.read_csv(path, nrows=0).columns
        sample_cols = cols[~cols.isin(index_cols)]

        ## Load feature values
        matrix = pd.read_csv(
            path,
            low_memory=False,
            usecols=sample_cols,
            dtype="float32",
            engine="c"
        )

        matrix = matrix.transpose()
        matrix.fillna(0, inplace=True)
        matrix = np.log1p(matrix)

        ## Assign feature names
        feature_names_cols = pd.read_csv(path, low_memory=False, usecols=index_cols, dtype=str)
        matrix.columns = \
            feature_names_cols["Chromosome"] + ";" + \
            feature_names_cols["PosStart"] + ";" + \
            feature_names_cols["PosEnd"]

        self._add_prefix_to_columns(matrix, "alt_sj")

        return matrix

    ## Combine ================================
    @cached_property
    def dna_features(self):
        features = dict(
            gen_pos = self.gen_pos_matrix,
            snv96 = self.snv96_matrix,
            event = self.event_matrix,
            sig = self.sig_matrix
        )

        for name, df in features.items():
            self._add_prefix_to_columns(df, prefix=name)

        features = pd.concat(features.values(), axis=1)
        features.fillna(0, inplace=True)

        features.index.name = "sample_id"

        return features

    @cached_property
    def rna_features(self):
        features = pd.concat([
            self.gene_exp_matrix,
            self.alt_sj_matrix
        ], axis=1)

        features.index.name = "sample_id"

        return features

    def load_dna_features(self):
        return CuppaFeatures(self.dna_features)

    def load_rna_features(self):
        return CuppaFeatures(self.rna_features)

    def load_features(self):

        df = pd.concat([
            self.dna_features,
            self.rna_features
        ], axis=1)

        return CuppaFeatures(df)
