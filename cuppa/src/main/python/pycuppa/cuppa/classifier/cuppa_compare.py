import numpy as np
import pandas as pd
from typing import Literal
from pprint import pformat
from functools import cached_property
from cuppa.classifier.cuppa_prediction import CuppaPredSummary
from cuppa.logger import LoggerMixin
from cuppa.misc.utils import as_categorical


class CuppaCompare(LoggerMixin):

    """
    Compares two CuppaPredSummary objects

    This class can be instantiated directly given two CuppaPredSummary objects. One is referred to as 'old' and the
    other as 'new'. Alternatively, prediction summaries stored as tsv files can also be loaded using the
    from_pred_summ_files() method.

    The following properties contain the comparison stats:

        prediction_comparison: Comparison of the top-1 prediction probabilities and classes for each sample.
        The 'correct_type' column shows whether the prediction in old only, new only, both or neither.

        correct_type_stats: Counts of each 'correct_type' per classifier

        performance_comparison: Comparison of recall and precision per class. This table also shows a comparison per
        class of the number of samples in total, predicted as that class, and correctly predicted.
    """

    def __init__(
        self,
        pred_summ_old: CuppaPredSummary,
        pred_summ_new: CuppaPredSummary,
    ):
        self.pred_summ_old = pred_summ_old
        self.pred_summ_new = pred_summ_new

        self._check_actual_class_column_exists()


    def __repr__(self) -> str:
        return pformat(vars(self))

    def _check_actual_class_column_exists(self):
        error_msg = "`actual_class` column is required but absent from `%s`"

        if not self.pred_summ_old._has_actual_class_column:
            self.logger.error(error_msg % "pred_summ_old")
            raise KeyError

        if not self.pred_summ_new._has_actual_class_column:
            self.logger.error(error_msg % "pred_summ_new")
            raise KeyError

    @staticmethod
    def from_pred_summ_files(path_old: str, path_new: str):
        return CuppaCompare(
            pred_summ_old = CuppaPredSummary.from_tsv(path_old),
            pred_summ_new = CuppaPredSummary.from_tsv(path_new)
        )

    @staticmethod
    def _align_rows(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        index_columns: str | list[str] | pd.Series,
        index_location: Literal["both", "columns", "index"] = "both",
        copy: bool = True
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        if copy:
            df1 = df1.copy()
            df2 = df2.copy()

        ## Get all rows, given `index_columns`
        index_union = pd.concat([
            df1[index_columns],
            df2[index_columns],
        ])

        ## Make indexes
        if len(pd.Series(index_columns))==1:
            index_union = index_union.unique()

            df1.index = df1[index_columns]
            df2.index = df2[index_columns]
        else:
            index_union = index_union.drop_duplicates()
            index_union = pd.MultiIndex.from_frame(index_union)

            df1.index = pd.MultiIndex.from_frame(df1[index_columns])
            df2.index = pd.MultiIndex.from_frame(df2[index_columns])

        ## Align indexes
        df1 = df1.reindex(index_union)
        df2 = df2.reindex(index_union)

        if index_location == "columns":
            df1.reset_index(inplace=True, drop=True)
            df2.reset_index(inplace=True, drop=True)

        if index_location == "index":
            df1.drop(index_columns, inplace=True, axis=1)
            df2.drop(index_columns, inplace=True, axis=1)

        return df1, df2

    @staticmethod
    def _move_index_to_columns(df: pd.DataFrame):
        df.index.names = pd.MultiIndex.from_arrays([
            ["info"] * len(df.index.names),
            df.index.names
        ])
        df = df.reset_index()
        return df

    @cached_property
    def prediction_comparison(self) -> pd.DataFrame:

        new = self.pred_summ_new.copy()
        old = self.pred_summ_old.copy()

        required_columns = ["sample_id", "clf_name", "actual_class", "pred_class_1", "pred_prob_1", "is_correct_pred"]
        new = new[required_columns]
        old = old[required_columns]

        new, old = self._align_rows(
            df1 = new,
            df2 = old,
            index_columns = ["sample_id", "clf_name"],
            index_location = "index",
            copy = False
        )

        def get_correct_type():
            _new = new["is_correct_pred"]
            _old = old["is_correct_pred"]

            correct_type = np.full(len(new), np.nan, dtype=object)

            correct_type[(_new == True)  & (_old == True) ] = "both"
            correct_type[(_new == True)  & (_old == False)] = "new_only"
            correct_type[(_new == False) & (_old == True) ] = "old_only"
            correct_type[(_new == False) & (_old == False)] = "neither"

            return pd.Series(correct_type, index=new.index)

        def column_equal(column: str) -> pd.Series:
            output = new[column] == old[column]
            output[new[column].isna() | old[column].isna()] = np.nan
            return output

        info = pd.DataFrame(dict(
            correct_type=get_correct_type()
        ))

        is_equal = pd.DataFrame(dict(
            actual_class = column_equal("actual_class"),
            pred_class_1 = column_equal("pred_class_1"),
            pred_prob_1 = column_equal("pred_prob_1"),
        ))

        comparison = pd.concat(
            dict(info=info, is_equal=is_equal, new=new, old=old),
            axis=1
        )

        comparison = self._move_index_to_columns(comparison)

        ## Force classifier order
        clf_names = comparison[("info","clf_name")]
        clf_names = as_categorical(clf_names)
        comparison[("info", "clf_name")] = clf_names

        return comparison

    @cached_property
    def correct_type_stats(self) -> pd.DataFrame:
        comparison = self.prediction_comparison

        stats = comparison["info"]\
            .groupby("clf_name")\
            ["correct_type"]\
            .value_counts(dropna=False)\
            .reset_index()

        return stats

    @cached_property
    def performance_comparison(self) -> pd.DataFrame:
        new = self.pred_summ_new.performance()
        old = self.pred_summ_old.performance()

        new, old = self._align_rows(
            df1=new,
            df2=old,
            index_columns=["class", "clf_name"],
            index_location="index",
            copy=False
        )

        diff = new - old

        comparison = pd.concat(
            dict(diff=diff, new=new, old=old),
            axis=1
        )

        comparison = self._move_index_to_columns(comparison)

        return comparison



