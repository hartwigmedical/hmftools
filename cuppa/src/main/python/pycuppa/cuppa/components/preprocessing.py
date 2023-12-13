import numpy as np
import pandas as pd
from typing import Self, Any, Optional
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted
from numpy.typing import NDArray

from cuppa.logger import LoggerMixin


class Log1pTransformer(BaseEstimator, TransformerMixin):
    def __init__(self):
        """
        ln(x+1) transformation
        """
        pass

    def fit(self, X: Any, y = None) -> Self:
        return self

    def transform(self, X: pd.DataFrame, y = None) -> pd.DataFrame | NDArray:
        return np.log1p(X)

    def set_output(self, transform = None) -> Self:
        return self


class NaRowFilter(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        pattern: Optional[str] = None,
        use_first_col: bool = False,
        show_warnings: bool = True
    ):
        """

        Parameters
        ----------
        pattern : str, default=None
            Name of columns containing this regex pattern will be included. If
            None, column selection will not be selected based on pattern.

        use_first_col : bool
            If true, assumes that NA rows found in the first matched column apply to all rows. This is done to speed up
            dataframe subsetting.

        show_warnings : bool
            Print number of NA rows that were removed?

        """

        self.pattern = pattern
        self.use_first_col = use_first_col
        self.show_warnings = show_warnings

    def fit(self, X: Any, y: None) -> Self:
        return self

    @staticmethod
    def detect_na_rows(
        X: pd.DataFrame,
        pattern: Optional[str] = None,
        use_first_col: bool = False,
    ) -> NDArray:

        ## Select columns to parse
        columns = X.columns

        if len(columns) > 0:
            if pattern is not None:
                columns = columns[columns.str.match(pattern)]

            if use_first_col:
                columns = columns[0]

        columns = pd.Series(columns)

        ## Find NAs
        is_na_row = np.isnan(X[columns]).any(axis=1)

        return is_na_row

    def transform(
        self,
        X: pd.DataFrame,
        y: Optional[pd.Series] = None
    ) -> pd.DataFrame | tuple[pd.DataFrame, pd.Series]:

        is_na_row = self.detect_na_rows(
            X=X,
            pattern=self.pattern,
            use_first_col=self.use_first_col
        )

        n_na_rows = np.sum(is_na_row)

        if n_na_rows == 0:
            if y is None:
                return X
            return X, y

        if self.show_warnings:
            self.logger.debug("Removed %s rows with NAs" % str(n_na_rows))

        if y is None:
            return X[~is_na_row]
        return X[~is_na_row], y[~is_na_row]

    def fit_transform(self, X, y: Optional[pd.Series] = None) -> pd.DataFrame:
        return self.transform(X, y)

    def set_output(self, transform = None) -> Self:
        return self


class MaxScaler(BaseEstimator, TransformerMixin, LoggerMixin):

    def __init__(self, clip: bool = True):
        """
        Scales values from 0 to 1 by dividing by the maximum value in the training set.


        Parameters
        ----------
        clip : bool
            Set values <0 and >1 after division to 0 and 1 respectively?

        """
        self.clip = clip

    def fit(self, X: pd.DataFrame, y: pd.Series = None) -> Self:

        max_values = np.max(X, axis=0)
        max_values = pd.Series(max_values, index=X.columns)
        self.max_values = max_values

        return self

    def transform(self, X: pd.DataFrame, y: pd.Series = None) -> pd.DataFrame | NDArray:

        check_is_fitted(self, "max_values")

        index = X.index
        columns = X.columns

        with np.errstate(divide='ignore', invalid='ignore'): ## Ignore divide by zero warnings
            X = np.divide(np.array(X), np.array(self.max_values))

        X = np.nan_to_num(X, nan=0) ## Convert NAs from divide by zero

        if self.clip:
            X = np.clip(X, a_min=0, a_max=1)

        X = pd.DataFrame(X, index=index, columns=columns)

        return X

    def set_output(self, transform = None) -> Self:
        return self





