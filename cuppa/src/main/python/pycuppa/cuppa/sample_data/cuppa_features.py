from __future__ import annotations

import os.path
from functools import cached_property
from typing import Iterable, Literal

import numpy as np
import pandas as pd

from cuppa.components.preprocessing import NaRowFilter
from cuppa.logger import LoggerMixin
from cuppa.constants import RESOURCES_DIR
from cuppa.misc.utils import check_required_columns


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
    ) -> "CuppaFeaturesPaths":

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

    def get_feat_type_cols(self, sel_feat_types: str | Iterable[str]) -> pd.Series:
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
    def from_tsv_files(cls, paths: CuppaFeaturesPaths, verbose: bool = True) -> "CuppaFeatures":

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

    def __init__(
        self,
        path: str,
        sample_id: str | None = None,
        excl_chroms: Iterable[str] | None = ["ChrY", "Y"],
        verbose: bool = True
    ):
        self.path = path
        self.sample_id = sample_id
        self.excl_chroms = excl_chroms
        self.verbose = verbose

    REQUIRED_COLS = ["Source","Category","Key","Value"]

    def load_file(self) -> pd.DataFrame:
        df = pd.read_table(self.path)
        check_required_columns(df, self.REQUIRED_COLS)

        if self.verbose:
            self.logger.info("Loaded file from: " + self.path)

        return df

    def add_sample_ids(self, df: pd.DataFrame) -> pd.DataFrame:

        if "SampleId" in df.columns:
            return df

        if self.sample_id is not None:
            df["SampleId"] = self.sample_id
            return df

        df["SampleId"] = "SAMPLE_1"
        self.logger.warning("No `sample_id` provided. Using `sample_id`='SAMPLE_1'")
        return df

    def filter_gen_pos_chroms(self, df: pd.DataFrame) -> pd.DataFrame:

        if self.excl_chroms is None:
            return df

        regex = "^(" + "|".join(self.excl_chroms) + ")"
        is_excluded_feature = (df["Category"] == "gen_pos") & df["Key"].str.match(regex)

        if self.verbose:
            df_excluded = df[is_excluded_feature]
            if len(df_excluded) > 0:
                self.logger.info(
                    "Removed %i gen_pos bin(s) in %i sample(s) with the regex '%s'",
                    len(df_excluded["Key"].unique()),
                    len(df_excluded["SampleId"].unique()),
                    regex
                )

        return df[~is_excluded_feature]

    def dedup_drivers(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.sort_values(by="Value", ascending=False)

        columns_orig = df.columns

        sample_feature_keys = pd.Series(df[["SampleId","Category","Key"]].itertuples(index=False))
        df.loc[:,"is_dup_feature"] = sample_feature_keys.duplicated().values
        df.loc[:,"is_driver"] = df["Category"] == "event.driver"

        df_dup = df.query("is_dup_feature & is_driver")
        if len(df_dup) > 0:
            self.logger.warning(
                "Removed %i duplicate driver(s) from %i sample(s) (kept highest impact driver)",
                len(df_dup), len(df_dup["SampleId"].unique())
            )

        is_excluded_feature = df["is_dup_feature"] & df["is_driver"]
        df = df[~is_excluded_feature]

        df = df.sort_index()
        return df[columns_orig]

    MAPPINGS_SIGS = dict(
        SIG_1="Age (SBS1)",
        SIG_2_13_AID_APOBEC="AID/APOBEC (SBS2/13)",
        SIG_4_SMOKING="Smoking (SBS4)",
        SIG_6_MMR="MMRD (SBS6)",
        SIG_7_UV="UV (SBS7)",
        SIG_10_POLE="POLE (SBS10)",
        SIG_11="Temozolomide (SBS11)",
        SIG_17="ROS/5FU (SBS17)",
    )

    EXCLUDED_SIGS = (
        "Age (SBS1)",
        "POLE (SBS10)",
        "Temozolomide (SBS11)"
    )

    def parse_signatures(self, df: pd.DataFrame) -> pd.DataFrame:
        df["Key"] = df["Key"].replace(self.MAPPINGS_SIGS)
        df = df[~df["Key"].isin(self.EXCLUDED_SIGS)]
        return df

    def reshape_long_to_wide(self, df: pd.DataFrame):

        df["FeatureNames"] = df["Category"] + "." + df["Key"]

        df_wide = df.pivot_table(index="SampleId", columns="FeatureNames", values="Value")
        df_wide = df_wide[df["FeatureNames"].unique()] ## Preserve original feature order

        df_wide.index.name = "sample_id"
        df_wide.columns.name = None

        df_wide = df_wide.fillna(0)

        return df_wide

    def load(self) -> pd.DataFrame:
        df = self.load_file()
        df = self.add_sample_ids(df)

        df = self.filter_gen_pos_chroms(df)
        df = self.dedup_drivers(df)
        df = self.parse_signatures(df)

        features = self.reshape_long_to_wide(df)

        n_samples = len(df["SampleId"].unique())
        n_features = features.shape[1]
        if self.verbose:
            if n_samples == 0:
                self.logger.info("Loaded %i features from sample_id: %s" % (n_features, df["SampleId"][0]))
            else:
                self.logger.info("Loaded %i features from %i samples" % (n_features, n_samples))

        return CuppaFeatures(features)


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

        feat_names = pd.read_csv(RESOURCES_DIR / ("feature_names/gen_pos.hg%i.csv" % self.genome_version))

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

        feat_names = pd.read_csv(RESOURCES_DIR / "feature_names/snv96.csv")
        matrix.columns = feat_names["context"].values

        return matrix

    @cached_property
    def sig_matrix(self) -> pd.DataFrame:

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
        df["sample_feature_id"] = df["SampleId"] + "_" + df["Name"] + "_" + df["sub_type"]
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
    def _df_trait(self) -> pd.DataFrame:
        return pd.read_csv(self.paths["trait"], index_col="SampleId")

    @cached_property
    def tmb_matrix(self) -> pd.DataFrame:
        matrix = pd.DataFrame({
            "snv_count": self.snv96_matrix.sum(axis=1),
            "indels_per_mb": self._df_trait["MsIndelsPerMb"]
        })
        self._add_prefix_to_columns(matrix, "tmb")
        return matrix

    @cached_property
    def trait_matrix(self) -> pd.DataFrame:
        matrix = pd.DataFrame({
            "is_male": (self._df_trait["Gender"] == "MALE").astype(int),
            "whole_genome_duplication": self._df_trait["WholeGenomeDuplication"].astype(int)
        })
        self._add_prefix_to_columns(matrix, "trait")
        return matrix


    @cached_property
    def event_matrix(self) -> pd.DataFrame:

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
    def gene_exp_matrix(self) -> pd.DataFrame:

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
    def alt_sj_matrix(self) -> pd.DataFrame:

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
    def dna_features(self) -> pd.DataFrame:
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
    def rna_features(self) -> pd.DataFrame:
        features = pd.concat([
            self.gene_exp_matrix,
            self.alt_sj_matrix
        ], axis=1)

        features.index.name = "sample_id"

        return features

    def load_dna_features(self) -> CuppaFeatures:
        return CuppaFeatures(self.dna_features)

    def load_rna_features(self) -> CuppaFeatures:
        return CuppaFeatures(self.rna_features)

    def load_features(self) -> CuppaFeatures:

        df = pd.concat([
            self.dna_features,
            self.rna_features
        ], axis=1)

        return CuppaFeatures(df)
