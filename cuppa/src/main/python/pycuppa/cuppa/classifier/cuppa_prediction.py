from __future__ import annotations

import logging
from functools import cached_property
from typing import Iterable, TYPE_CHECKING

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pandas.core.dtypes.common import is_integer_dtype
from sklearn.compose import make_column_selector

from cuppa.constants import CLF_GROUPS, SUB_CLF_NAMES
from cuppa.logger import LoggerMixin
from cuppa.misc.utils import get_top_cols, check_required_columns, as_categorical
from cuppa.performance.performance_stats import PerformanceStatsBuilder, PerformanceStats

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from cuppa.classifier.cuppa_classifier import CuppaClassifier


class CuppaPredictionBuilder(LoggerMixin):

    def __init__(
        self,
        cuppa_classifier: CuppaClassifier,
        X: pd.DataFrame,
        verbose: bool = False
    ):
        self.cuppa_classifier = cuppa_classifier
        self.X = X
        self.verbose = verbose

    INDEX_NAMES = ["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value"]

    @property
    def class_columns(self) -> list[str]:
        return list(self.cuppa_classifier.classes_)

    @property
    def all_columns(self) -> list[str]:
        return self.INDEX_NAMES + self.class_columns

    def _set_required_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.reset_index().reindex(self.all_columns, axis=1)

    def _set_required_indexes(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.set_index(self.INDEX_NAMES)

    @cached_property
    def probs(self) -> pd.DataFrame:
        if self.verbose:
            self.logger.info("Getting probabilities")

        df = self.cuppa_classifier.predict_proba(self.X)
        df["data_type"] = "prob"

        df = self._set_required_columns(df)
        df = self._set_required_indexes(df)

        return df

    @cached_property
    def feat_contribs(self) -> pd.DataFrame:
        if self.verbose:
            self.logger.info("Getting feature contributions")

        ## Contribs --------------------------------
        df = self.cuppa_classifier.feat_contrib(X=self.X, sub_clf_names=SUB_CLF_NAMES.EVENT, column_type="classes")
        df = self._set_required_columns(df)

        ## Row info
        df["data_type"] = "feat_contrib"
        df["clf_group"] = CLF_GROUPS.DNA

        ## Get feat values --------------------------------
        feat_values = self.X[make_column_selector("^event")]

        ## Make single column dataframe
        feat_values = feat_values.stack().to_frame(name="feat_value")
        feat_values.index.names = ("sample_id", "feat_name")

        ## Align `feat_values.index` with `feat_contribs.index
        index = pd.MultiIndex.from_frame(df[["sample_id","feat_name"]])
        feat_values = feat_values.reindex(index)

        df["feat_value"] = feat_values.values

        df = self._set_required_indexes(df)

        return df

    @cached_property
    def sig_quantiles(self) -> pd.DataFrame:
        if self.verbose:
            self.logger.info("Getting mutational signature quantiles")

        df = self.cuppa_classifier.sig_quantile_transformer.transform(
            self.X[make_column_selector("^sig")], include_feat_values=True
        )

        df = self._set_required_columns(df)

        ## Add row info
        df["data_type"] = "sig_quantile"
        df["feat_name"] = df["feat_name"].replace("^\w+[.]", "", regex=True)

        df = self._set_required_indexes(df)

        return df


    def build(self) -> CuppaPrediction:

        df = pd.concat([self.probs, self.feat_contribs, self.sig_quantiles])
        df = df.loc[self.X.index] ## Make rows from the same sample consecutive

        return CuppaPrediction.from_data_frame(df)


class CuppaPrediction(pd.DataFrame, LoggerMixin):
    def __init__(self, predictions: pd.DataFrame, *args, **kwargs):
        super().__init__(predictions, *args, **kwargs)

    @property
    def _constructor(self):
        return CuppaPrediction

    ## ================================
    @property
    def class_columns(self) -> list[str]:
        return list(self.columns)

    @cached_property
    def sample_ids(self) -> pd.Index:
        return self.index.get_level_values("sample_id").dropna().unique()

    @property
    def is_multi_sample(self) -> bool:
        return len(self.sample_ids.dropna()) > 1

    ## Utility methods ================================
    @staticmethod
    def concat(predictions_list) -> CuppaPrediction:
        predictions_merged = pd.concat(predictions_list, axis=0)
        return CuppaPrediction(predictions_merged)

    ## Subtypes / supertypes ================================
    @staticmethod
    def _columns_to_class_metadata(columns: Iterable[str]) -> pd.DataFrame:

        classes = pd.Series(columns)

        ## Get subtypes / supertypes
        prefixes = classes.str.extract("(^[A-Za-z/ ]+)", expand=False)
        suffixes = classes.str.replace("^[A-Za-z/ ]+: ", "", regex=True)

        metadata = pd.DataFrame(
            dict(
                supertype = prefixes.values,
                subtype = suffixes.values,
            ),
            index=classes.values
        )

        ## Sort cancer types, with "Other" subtypes at the end
        metadata["is_other_subtype"] = metadata["subtype"].str.endswith("Other")
        metadata = metadata.sort_values(["supertype", "is_other_subtype", "subtype"])

        ## Count number of subtypes in supertype
        supertype_counts = metadata["supertype"].value_counts()
        metadata["subtypes_in_supertype"] = supertype_counts[metadata["supertype"]].values

        return metadata

    @cached_property
    def class_metadata(self) -> pd.DataFrame:
        return self._columns_to_class_metadata(self.columns)

    @property
    def supertypes(self) -> pd.Series:
        return self.class_metadata["supertype"]


    ## I/O ================================
    @classmethod
    def from_data_frame(cls, df: pd.DataFrame) -> CuppaPrediction:
        check_required_columns(df, required_columns=CuppaPredictionBuilder.INDEX_NAMES, from_index=True)

        ## Force column order with Other subtypes as last amongst subtypes
        metadata = cls._columns_to_class_metadata(df.columns)
        df = df[metadata.index]

        return CuppaPrediction(df)

    def to_tsv(
        self,
        path: str,
        signif_digits: int | None = 4,
        chunksize: int = 1000,
        verbose: bool = False,
        *args, **kwargs
    ) -> None:

        ## Use two row column headers to mark index columns and supertype columns
        columns = list(self.index.names) + list(self.columns)
        column_metadata = ["index"]*len(self.index.names) + list(self.supertypes.values)

        df = self.reset_index()
        df.columns = pd.MultiIndex.from_arrays([column_metadata, columns])

        float_format = f"%.{signif_digits}g" if signif_digits is not None else None

        if verbose:
            self.logger.info("Writing predictions to: " + path)

        df.to_csv(
            path,
            sep="\t",
            index=False,
            float_format=float_format,
            chunksize=chunksize,
            *args, **kwargs
        )

    @classmethod
    def from_tsv(cls, path: str) -> CuppaPrediction:
        df = pd.read_csv(path, sep="\t", header=[0,1])

        index_cols = df["index"].columns
        df.columns = df.columns.get_level_values(1)
        df.set_index(list(index_cols), inplace=True)

        return cls.from_data_frame(df)

    def subset_probs_by_clf_groups(self, clf_groups: str | list[str] = None, verbose: bool = True) -> CuppaPrediction:

        clf_groups = pd.Series(clf_groups)
        if not all(clf_groups.isin(CLF_GROUPS.get_all())):
            self.logger.error("`clf_groups` must be one or more of: " + ", ".join(CLF_GROUPS.get_all()))
            raise ValueError

        if verbose:
            self.logger.info("Subsetting probabilties for clf_groups: " + ", ".join(clf_groups))

        indexes = self.index.to_frame(index=False)
        selected_rows = (
            (indexes["data_type"] != "prob") |
            (
                (indexes["data_type"] == "prob") &
                (indexes["clf_group"].isin(clf_groups))
            )
        )

        return self[selected_rows.values]

    def wide_to_long(self) -> pd.DataFrame:

        self.columns.name = "cancer_type"
        series = self.stack(dropna=False).to_frame("data_value")
        self.columns.name = None

        long = series.index.to_frame(index=False)

        long["data_value"] = series.values

        return long

    ## Get specific data types ================================
    def get_rows(self, keys: str | Iterable[str], index_level: str = "sample_id") -> CuppaPrediction:
        keys = pd.Series(keys)

        index_level_values = self.index.get_level_values(index_level)
        is_selected_row = index_level_values.isin(keys)

        if not is_selected_row.any():
            self.logger.error("Some `keys` were not found in index level `%s`" % index_level)
            raise LookupError

        return self.loc[is_selected_row]

    @cached_property
    def probs(self) -> CuppaPrediction:
        probs = self.get_rows("prob", index_level="data_type")
        probs.index = probs.index.droplevel(["data_type", "feat_name", "feat_value"])
        return probs

    def get_data_types(self, keys: str | Iterable[str]) -> CuppaPrediction:
        return self.get_rows(keys, index_level="data_type")

    def get_clf_names(self, keys: str | Iterable[str]) -> CuppaPrediction:
        return self.get_rows(keys, index_level="clf_name")

    def get_clf_groups(self, keys: str | Iterable[str]) -> CuppaPrediction:
        return self.get_rows(keys, index_level="clf_group")

    def get_samples(self, keys: str | int | Iterable[str | int]) -> CuppaPrediction:
        #keys=[0,2]

        keys = pd.Series(keys)
        if is_integer_dtype(keys):
            keys = self.sample_ids[keys]

        return self.get_rows(keys=keys, index_level="sample_id")

    ## Downstream output ================================
    def summarize(
        self,
        actual_classes: pd.Series | None = None,
        top_n_classes: int = 3,
        top_n_features: int = 3,
        show_extra_info: bool = False,
        verbose: bool = False
    ) -> CuppaPredSummary:

        builder = CuppaPredSummaryBuilder(
            predictions = self,
            actual_classes = actual_classes,

            top_n_features = top_n_features,
            top_n_classes = top_n_classes,
            show_extra_info= show_extra_info,
            verbose = verbose
        )

        pred_summ = builder.build()

        return pred_summ


class CuppaPredSummaryBuilder(LoggerMixin):
    def __init__(
        self,
        predictions: CuppaPrediction,
        actual_classes: pd.Series | None = None,

        top_n_features: int = 3,
        top_n_classes: int = 3,
        show_extra_info: bool = False,
        verbose: bool = False
    ):
        self.predictions = predictions
        self.actual_classes = actual_classes

        self.top_n_features = top_n_features
        self.top_n_classes = top_n_classes
        self.show_extra_info = show_extra_info
        self.verbose = verbose

    ## Top predictions ================================
    @cached_property
    def _sorted_preds(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        sorted_classes, sorted_probs = get_top_cols(self.predictions.probs)
        return sorted_classes, sorted_probs

    @cached_property
    def _top_preds(self) -> tuple[pd.DataFrame, pd.DataFrame]:

        if self.verbose:
            self.logger.debug("Getting top cancer types")

        sorted_classes, sorted_probs = self._sorted_preds
        top_classes = sorted_classes.iloc[:, 0:self.top_n_classes]
        top_probs = sorted_probs.iloc[:, 0:self.top_n_classes]
        return top_classes, top_probs

    @property
    def top_classes(self) -> pd.DataFrame:
        return self._top_preds[0]

    @property
    def top_probs(self) -> pd.DataFrame:
        return self._top_preds[1]

    ## Top signatures ================================
    @staticmethod
    def _collapse_columns_to_strings(df: pd.DataFrame, sep=", ") -> pd.Series:

        strings = df.iloc[:,0].copy()
        for i in range(len(df.columns)):
            if i == 0: continue
            column = df.iloc[:,i].copy()
            column_with_sep = sep + column
            column_with_sep[column.str.len() == 0] = ""
            strings += column_with_sep

        return strings


    @cached_property
    def tmb_snv_per_sample(self) -> pd.Series:
        index = self.predictions.get_data_types("feat_contrib").index.to_frame(index=False)
        index = index[index["feat_name"].str.endswith("snv_count")]

        tmb = index["feat_value"].astype(int)
        tmb.index = index["sample_id"]
        return tmb

    @cached_property
    def extra_info_sigs(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting top signatures")

        ## Get signature info
        index = self.predictions \
            .get_data_types("sig_quantile") \
            .index.to_frame(index=False)

        ## Get top signatures
        index["rank"] = index.groupby("sample_id")["feat_value"] \
            .rank(method="first", ascending=False) \
            .astype(int)

        index["sample_id"] = as_categorical(index["sample_id"])
        index = index.sort_values(["sample_id", "rank"])

        ## Get relative contribution
        index["tmb"] = self.tmb_snv_per_sample[index["sample_id"]].values
        index["perc"] = (index["feat_value"] / index["tmb"]) * 100

        index["clf_name"] = "snv96"
        index["clf_group"] = "snv96"

        ## Make strings
        sig_names = index.pivot(index="sample_id", columns="rank", values="feat_name")
        sig_values = index.pivot(index="sample_id", columns="rank", values="feat_value").round(1)
        sig_percs = index.pivot(index="sample_id", columns="rank", values="perc").round(1)

        labels = sig_names + "=" + sig_values.astype(str) + "," + sig_percs.astype(str) + "%"
        labels[sig_values == 0] = ""

        ##
        strings = self._collapse_columns_to_strings(labels, "; ")
        is_empty_string = strings.str.len() == 0
        format = np.where(is_empty_string, "", "{sig_name}={count},{percent}")

        df = pd.DataFrame(dict(
            extra_info=strings,
            extra_info_format=format
        ))

        df.index = pd.MultiIndex.from_frame(pd.DataFrame(dict(
            sample_id=sig_names.index,
            clf_group="dna",
            clf_name="snv96"
        )))

        return df

    ## Top features ================================
    @cached_property
    def _top_features(self) -> tuple[pd.DataFrame, pd.DataFrame]:

        ## Get event features
        if self.verbose:
            self.logger.debug("Getting feature contributions")

        contribs = self.predictions.get_rows("feat_contrib", index_level="data_type")

        ## Simplify index
        contribs.index = pd.MultiIndex.from_arrays([
            contribs.index.get_level_values("sample_id"),
            contribs.index.get_level_values("feat_name")
        ])

        ## Drop prior
        contribs = contribs.loc[contribs.index.get_level_values("feat_name") != "_prior"]

        ## Transpose: (sample_id, feat_name) x pred_class --> (sample_id, pred_class) x feat_name
        contribs.columns.name = "pred_class"
        contribs_T = contribs.stack().unstack("feat_name")

        ## Subset for contributions of top classes
        top_classes_long = self.top_classes.stack()\
            .to_frame("pred_class")\
            .reset_index()\
            .query("clf_name=='event'")

        selected_indexes = pd.MultiIndex.from_frame(top_classes_long[["sample_id", "pred_class"]])
        contribs_T = contribs_T.loc[selected_indexes]

        ## Get top features
        top_features, top_contribs = get_top_cols(contribs_T, top_n=self.top_n_features)

        return top_features, top_contribs

    @property
    def top_features(self) -> pd.DataFrame:
        return self._top_features[0]

    @property
    def top_contribs(self) -> pd.DataFrame:
        return self._top_features[1]

    @staticmethod
    def _signif(df, precision=2) -> pd.DataFrame:
        index = df.index
        columns = df.columns

        arr = np.asarray(df)

        x_positive = np.where(np.isfinite(arr) & (arr != 0), np.abs(arr), 10 ** (precision - 1))
        mags = 10 ** (precision - 1 - np.floor(np.log10(x_positive)))
        arr = np.round(arr * mags) / mags

        return pd.DataFrame(arr, index=index, columns=columns)

    @cached_property
    def top_feat_values(self) -> pd.DataFrame:
        ## Extract feature values from index of predictions
        feat_values = self.predictions.index.to_frame(index=False)

        ## Subset for relevant columns
        feat_values = feat_values[["sample_id", "feat_name", "feat_value"]]

        ## Set index
        feat_values = feat_values.set_index(["sample_id", "feat_name"])

        ## Align indexes of
        top_features_long = self.top_features.stack().to_frame("feat_name").reset_index()

        selected_indexes = pd.MultiIndex.from_frame(top_features_long[["sample_id", "feat_name"]])
        feat_values = feat_values.loc[selected_indexes]

        ## Long to wide; results in same dimensions as `top_features`
        feat_values["rank"] = top_features_long["rank"].values
        feat_values["pred_class"] = top_features_long["pred_class"].values

        top_feat_values = feat_values \
            .reset_index() \
            .pivot(index=["sample_id", "pred_class"], columns="rank", values="feat_value")

        ## Align sample_id and pred_class
        top_feat_values = top_feat_values.loc[self.top_features.index]

        return top_feat_values

    @cached_property
    def extra_info_feat(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Getting top feature contributions and values")

        ## Dataframe: (sample_id, pred_class) x feat_rank
        strings_semi_long = \
            self.top_features + "=" + \
            self._signif(self.top_feat_values, precision=2).astype(str) + "," + \
            self.top_contribs.round(decimals=1).astype(str)

        ## Remove feat group prefix
        strings_semi_long = strings_semi_long.replace("^event[.]", "", regex=True)

        ## Series: (sample_id, pred_class)
        strings_long = self._collapse_columns_to_strings(strings_semi_long, sep='; ')

        ## Long to wide; result in sample dimensions as `top_classes`
        ## Dataframe: (sample_id, pred_class) x cancer_type_rank
        strings_long = strings_long.to_frame("string").reset_index()
        strings_long["rank"] = strings_long.groupby("sample_id").cumcount() + 1  ## Add rank
        strings_long["sample_id"] = as_categorical(strings_long["sample_id"])

        strings_wide = strings_long.pivot(index="sample_id", columns="rank", values="string")
        strings_wide = pd.Series(strings_wide.columns).astype(str).values + ": " + strings_wide

        ##
        df = pd.DataFrame(dict(
            extra_info=self._collapse_columns_to_strings(strings_wide, sep=' || '),
            extra_info_format="{pred_class_n}: {feat_name}={value},{log_odds}"
        ))

        df.index = pd.MultiIndex.from_frame(pd.DataFrame(dict(
            sample_id=strings_wide.index,
            clf_group="dna",
            clf_name="event"
        )))

        return df

    @property
    def extra_info(self) -> pd.DataFrame:
        return pd.concat([
            self.extra_info_sigs,
            self.extra_info_feat
        ])

    ## Compare prediction and actual class ================================
    @cached_property
    def correct_info(self) -> pd.DataFrame:

        if self.verbose:
            self.logger.debug("Determining whether each prediction is correct")

        ## Get actual class per row --------------------------------
        sample_ids = self.predictions.probs.index.get_level_values("sample_id")
        actual_classes = pd.Series(self.actual_classes)[sample_ids]

        ## Top prediction is correct? --------------------------------
        sorted_classes, sorted_probs = self._sorted_preds
        pred_classes = sorted_classes[1]

        is_correct_pred = actual_classes == pred_classes.values
        is_correct_pred[pred_classes.isna().values] = np.nan

        ## Which prediction is correct? --------------------------------
        ## Boolean matrix
        which_correct_pred = np.reshape(actual_classes.values, (-1, 1)) == sorted_classes.values

        ## Which column is true
        which_correct_pred = np.argmax(which_correct_pred, axis=1) + 1

        ## Propagate NAs
        which_correct_pred = pd.Series(which_correct_pred)
        which_correct_pred[pred_classes.isna().values] = np.nan

        correct_info = pd.DataFrame(
            dict(
                actual_class = actual_classes.values,
                is_correct_pred = is_correct_pred.values,
                which_correct_pred = which_correct_pred.values
            ),
            index=self.predictions.probs.index
        )

        # ## Left align sample_id and actual_class
        # correct_info = correct_info.reset_index().set_index(["sample_id", "actual_class"]).reset_index()

        return correct_info

    def build(self) -> CuppaPredSummary:

        if self.verbose:
            self.logger.info("Summarizing predictions")

        ## Main --------------------------------
        tables = dict()

        tables["pred_class"] = self.top_classes
        tables["pred_prob"] = self.top_probs

        ## Merge
        pred_summ = pd.concat(tables, axis=1)
        pred_summ.columns.names = ["result_type", "rank"]

        ## Flatten column multi index
        pred_summ.columns = \
            pred_summ.columns.get_level_values("result_type") + \
            "_" + \
            pred_summ.columns.get_level_values("rank").astype(str)

        ## Optional columns --------------------------------
        if self.show_extra_info:
            pred_summ = pd.concat([pred_summ, self.extra_info], axis=1)

        if self.actual_classes is not None:
            if self.verbose:
                self.logger.debug("Adding 'correct prediction' info")

            pred_summ = pd.concat([self.correct_info, pred_summ], axis=1)

        pred_summ.reset_index(inplace=True)

        if self.actual_classes is not None:
            ## Move actual_classes next to sample_id
            column_to_move = pred_summ.pop("actual_class")
            pred_summ.insert(1, "actual_class", column_to_move)

        return CuppaPredSummary.from_data_frame(pred_summ)


class CuppaPredSummary(pd.DataFrame, LoggerMixin):
    def __init__(self, pred_summ: pd.DataFrame, *args, **kwargs):
        super().__init__(pred_summ, *args, **kwargs)

    @property
    def _constructor(self):
        return CuppaPredSummary

    ## Init attributes ================================
    @cached_property
    def _has_actual_class_column(self) -> bool:
        return "actual_class" in self.columns

    @cached_property
    def classes(self) -> NDArray:
        classes = []

        classes += list(self["pred_class_1"].unique())
        if self._has_actual_class_column:
            classes += list(self["actual_class"].unique())

        classes = pd.Series(classes).sort_values().dropna().unique()
        return classes

    @cached_property
    def clf_names(self) -> NDArray:
        return self["clf_name"].unique()

    ## I/O ================================
    def to_tsv(self, path: str, verbose: bool = False) -> None:
        if verbose:
            self.logger.info("Writing prediction summary to: " + path)

        self.to_csv(path, sep="\t", index=False)

    @staticmethod
    def from_data_frame(df: pd.DataFrame) -> CuppaPredSummary:
        required_columns = ["sample_id", "clf_group", "clf_name", "pred_class_1", "pred_prob_1"]
        check_required_columns(df, required_columns=required_columns)

        ## Force categories
        df["sample_id"] = as_categorical(df["sample_id"])
        df["clf_name"] = as_categorical(df["clf_name"])
        df["clf_group"] = as_categorical(df["clf_group"])

        return CuppaPredSummary(df)

    @classmethod
    def from_tsv(cls, path: str) -> CuppaPredSummary:
        df = pd.read_table(path)
        return cls.from_data_frame(df)

    def long_to_wide(self) -> pd.DataFrame:

        ## Set constant/shared columns as index
        index_columns = ["sample_id"]
        if self._has_actual_class_column:
            index_columns += ["actual_class"]

        df = self.set_index(index_columns)

        ## Split by clf_name and concatenate horizontally
        df_split = {}
        for clf_name in self.clf_names:
            df_i = df[df["clf_name"] == clf_name]
            df_i = df_i.drop(["clf_group", "clf_name"], axis=1)  ## Remove and store classifier info columns
            df_split[clf_name] = df_i

        df = pd.concat(df_split, axis=1)

        ## Annotate column type for index
        df.index.names = pd.MultiIndex.from_arrays([
            ["metadata"] * len(df.index.names),
            df.index.names
        ])
        df = df.reset_index()

        return df

    def performance(self) -> PerformanceStats:
        return PerformanceStatsBuilder(self).build()

    def performance_by_prob_bin(
        self,
        bins: NDArray = np.linspace(0, 1, 6)
    ) -> PerformanceStats:
        return PerformanceStatsBuilder(self, prob_bins=bins).build(by_prob_bin=True)

