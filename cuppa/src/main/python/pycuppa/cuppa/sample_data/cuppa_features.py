from __future__ import annotations

import os.path
import pandas as pd

from cuppa.constants import NA_FILL_VALUE
from cuppa.logger import LoggerMixin
from cuppa.misc.utils import check_required_columns


class CuppaFeaturesDir(LoggerMixin):

    def __init__(self, directory: str):
        self.directory = directory

    FILE_PATTERNS = dict(
        snv="cuppa_data.cohort.snv.*.tsv",
        sv="cuppa_data.cohort.sv.*.tsv",
        trait="cuppa_data.cohort.sample_trait.*.tsv",
        driver="cuppa_data.cohort.driver.*.tsv",

        gene_exp="cuppa_data.cohort.gene_exp.*.tsv",
        alt_sj="cuppa_data.cohort.alt_sj.*.tsv",
    )

    OPTIONAL_PATTERN_KEYS = ("gene_exp", "alt_sj")

    @staticmethod
    def _find_files_in_dir_by_pattern(directory: str, pattern: str) -> pd.Series:
        files = pd.Series(os.listdir(directory))
        matched_files = files[files.str.match(pattern)]
        return matched_files

    def get_file_paths(self) -> pd.Series:

        paths = {}
        missing_any_required_files = False

        for key, pattern in self.FILE_PATTERNS.items():

            matched_files = self._find_files_in_dir_by_pattern(self.directory, pattern)

            if len(matched_files) == 0:

                if key in self.OPTIONAL_PATTERN_KEYS:
                    self.logger.warning("Missing optional input file type '%s' with pattern '%s'" % (key, pattern))
                else:
                    self.logger.error("Missing required input file type '%s' with pattern '%s'" % (key, pattern))
                    missing_any_required_files = True

                continue

            if len(matched_files) > 1:
                self.logger.warning("Pattern '%s' matched multiple files. Using the first: %s" % (pattern, ', '.join(matched_files)))

            matched_file = os.path.join(self.directory, matched_files.iloc[0])
            paths[key] = matched_file

        if missing_any_required_files:
            raise FileNotFoundError

        return pd.Series(paths)


class CuppaFeaturesLoader(LoggerMixin):
    def __init__(
        self,
        path: str,
        sample_id: str | None = None,
        na_fill_value: int | float = NA_FILL_VALUE
    ):
        self.path = path
        self.sample_id = sample_id
        self.na_fill_value = na_fill_value

        self.df: pd.DataFrame = None

    FLD_SOURCE = "Source"
    FLD_CATEGORY = "Category"
    FLD_KEY = "Key"
    FLD_VALUE = "Value"

    INDEX_COLS = [FLD_SOURCE, FLD_CATEGORY, FLD_KEY]

    CHROMS_AUTOSOMAL_AND_X = pd.Series([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]).astype(str)

    @property
    def is_multi_sample(self) -> bool:
        if os.path.isdir(self.path):
            return True
        else:
            header = pd.read_table(self.path, nrows=0).columns
            return len(header) > len(self.INDEX_COLS)+1

    @property
    def is_single_sample(self):
        return not self.is_multi_sample

    def _load_one_file(self) -> pd.DataFrame:
        df = pd.read_table(self.path)
        check_required_columns(df, self.INDEX_COLS)

        if self.is_single_sample:
            if self.sample_id is None:
                self.logger.error("`sample_id` must be provided when features are for one sample")
                raise ValueError

            df.rename(columns={self.FLD_VALUE: self.sample_id}, inplace=True)

        return df

    def _load_files_from_dir(self) -> pd.DataFrame:
        paths = CuppaFeaturesDir(directory=self.path).get_file_paths()
        dfs = {}
        for data_type, path in paths.items():
            self.logger.debug("Loading file: " + os.path.basename(path))
            df = pd.read_table(path)
            check_required_columns(df, self.INDEX_COLS)
            dfs[data_type] = pd.read_table(path)

        df_merged = pd.concat(dfs.values())

        return df_merged

    def _load_data(self) -> None:

        if os.path.isfile(self.path):

            self.logger.info(
                "Loading %s sample data file: %s",
                "single" if self.is_single_sample else "multi",
                self.path
            )

            self.df = self._load_one_file()

        else:
            self.logger.info("Loading multi sample data from dir: " + self.path)
            self.df = self._load_files_from_dir()

    def _assert_no_duplicate_features(self) -> None:
        duplicated_keys = self.df[self.FLD_KEY][self.df[self.FLD_KEY].duplicated()]
        if len(duplicated_keys):
            self.logger.error("The following features are duplicated: ", ", ".join(duplicated_keys.unique()))
            raise AssertionError

    def _check_gen_pos_chroms(self) -> None:
        bin_names = self.df.loc[ self.df[self.FLD_CATEGORY] == "gen_pos", "Key" ]
        chroms_uniq = pd.Series(bin_names.str.split("_").str[0].unique())

        unexpected_chroms = chroms_uniq[
            ~chroms_uniq.isin("chr" + self.CHROMS_AUTOSOMAL_AND_X) &
            ~chroms_uniq.isin(self.CHROMS_AUTOSOMAL_AND_X)
        ]

        if len(unexpected_chroms) > 0:
            self.logger.error("gen_pos bins must contain only chromosomes 1-22 and X. Found invalid chromosomes: " + ", ".join(unexpected_chroms))
            raise ValueError

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

    EXCLUDED_SIGS = ("SIG_1", "SIG_10_POLE", "SIG_11") ## TODO: Move this logic to java

    def _parse_signatures(self) -> None:
        self.df = self.df[~self.df[self.FLD_KEY].isin(self.EXCLUDED_SIGS)]
        self.df[self.FLD_KEY].replace(self.MAPPINGS_SIGS, inplace=True)

    def _assign_feature_names(self) -> None:
        feature_names = self.df[self.FLD_CATEGORY] + "." + self.df[self.FLD_KEY]
        self.df.drop(columns=self.INDEX_COLS, inplace=True)
        self.df.index = feature_names

    def _print_stats(self) -> None:

        n_samples = self.df.shape[1]
        n_features = self.df.shape[0]

        if self.is_multi_sample:
            self.logger.debug(f"Loaded {n_features} features from {n_samples} samples")
        else:
            self.logger.debug(f"Loaded {n_features} features from sample({self.sample_id})")

    def load(self) -> pd.DataFrame:
        self._load_data()
        self._assert_no_duplicate_features()
        self._check_gen_pos_chroms()
        self._parse_signatures()
        self._assign_feature_names()
        self._print_stats()

        return self.df.fillna(self.na_fill_value).transpose()
