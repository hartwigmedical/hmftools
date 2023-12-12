import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted
from typing import Self, Any, Literal
from numpy.typing import NDArray
from cuppa.logger import LoggerMixin
import logging
logger = logging.getLogger(__name__)

## Cosine similarity ================================
def cos_sim(
    A: NDArray | pd.Series | pd.DataFrame,
    B: NDArray | pd.Series | pd.DataFrame,
    return_pandas: bool = True
) -> NDArray | pd.DataFrame:

    """

    Parameters
    ----------
    A : array of shape (n_features) or matrix of shape (*, n_features)

    B : array of shape (n_features) or matrix of shape (*, n_features)
        If `B` is an array, `B` must have the same length as `A`
        If `B` is a matrix, `B` must have the same number of columns as `A`

    return_pandas : bool

    Returns
    -------

    """

    if A.ndim == 2:
        if B.ndim == 2:
            return _cos_sim_matrix(A, B, return_pandas=return_pandas)
        if B.ndim == 1:
            return _cos_sim_matrix(A, np.array(B).reshape(1, -1), return_pandas=return_pandas)

    if A.ndim == 1:
        if B.ndim == 1:
            return _cos_sim_array(A, B)
        if B.ndim == 2:
            return _cos_sim_matrix(np.array(A).reshape(1, -1), B, return_pandas=return_pandas)

    logger.error("A and B must be 1 or 2 dimensional arrays")
    raise TypeError()

def _cos_sim_array(a: NDArray | pd.Series, b: NDArray | pd.Series) -> NDArray:
    numerator = np.dot(a, b)
    denominator = np.sqrt(np.sum(a ** 2)) * np.sqrt(np.sum(b ** 2))

    with np.errstate(invalid='ignore'):
        result = numerator / denominator

    if np.isnan(result):
        return 0

    return result

def _cos_sim_matrix(
    A: NDArray | pd.Series | pd.DataFrame,
    B: NDArray | pd.Series | pd.DataFrame,
    return_pandas: bool = True
) -> NDArray | pd.DataFrame:

    numerator = np.dot(A, B.transpose())

    A_norm_L2 = np.sqrt(np.sum(A**2, axis=1))
    B_norm_L2 = np.sqrt(np.sum(B**2, axis=1))
    denominator = np.outer(A_norm_L2, B_norm_L2)

    with np.errstate(invalid='ignore'):
        M = numerator / denominator
    M[np.isnan(M)] = 0

    if return_pandas:
        M = pd.DataFrame(M)
        if isinstance(A, pd.DataFrame): M.index = A.index
        if isinstance(B, pd.DataFrame): M.columns = B.index

    return M


class ProfileSimilarityTransformer(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        agg_func: Literal["sum","mean","median","iqm"] = "sum",
        count_ceiling: int | float = None,
        normalize_profiles: bool = False,
        feature_prefix: str = None
    ):
        """
        Calculate similarity to per-cancer type mutation profiles

        Parameters
        ----------
        agg_func: str, 'sum', 'mean', 'iqm' (inter-quartile mean), or 'median'
            How to aggregate mutation profiles of the samples in each cancer type

        count_ceiling: int
            Maximum total mutation count per sample. If a sample has more mutations than `count_ceiling`, the array
            will be scaled such that the sum of the array equals `count_ceiling`. This is so that hyper-mutator samples
            (i.e. those with extremely high number of mutations) do not dominate the cancer type profile.

        normalize_profiles: bool
            If True, the resulting per-cancer type profiles will be normalized from 0 to 1

        metric: callable, default: cos_sim
            A function used to calculate the similarity of a sample to each per-cancer type profile

        feature_prefix: str
            A prefix to append to the per cancer type profile matrix


        """

        self.agg_func = agg_func
        self.count_ceiling = count_ceiling
        self.normalize_profiles = normalize_profiles
        self.feature_prefix = feature_prefix

    def _check_X_y_shape(self, X: pd.DataFrame | NDArray, y = pd.Series | NDArray):
        if X.ndim != 2:
            self.logger.error("X must be a 2D array-like")
            raise TypeError

        if y.ndim != 1:
            self.logger.error("y must be a 1D array-like")
            raise TypeError

        if X.shape[0] != y.shape[0]:
            self.logger.error("No. of X rows must equal length of y")
            raise Exception

    def _calc_column_averages(self, X):
        if self.agg_func == "sum":
            return np.sum(X, axis=0)

        if self.agg_func == "mean":
            return np.mean(X, axis=0)

        if self.agg_func == "median":
            return np.median(X, axis=0)

        if self.agg_func == "iqm":
            return trim_mean(X, proportiontocut=0.25, axis=0)

        self.logger.error("Invalid `agg_func`")
        raise ValueError()

    def _apply_count_ceiling(self, X):
        row_sums = X.sum(axis=1)
        row_ceilings = np.where(row_sums > self.count_ceiling, self.count_ceiling, row_sums)
        row_scale_factors = row_ceilings / row_sums

        X = X.multiply(row_scale_factors, axis=0)

        return X

    def fit(self, X: pd.DataFrame, y: pd.Series) -> Self:
        """

        Parameters
        ----------
        X: pandas DataFrame of shape (n_samples, n_features)
            Where `n_samples` is the number of samples and `n_features` is the number of features

        y: array-like of shape (n_samples,) default=None
            Labels for each sample

        Returns
        -------
        self
            Fitted transformer

        """

        self._check_X_y_shape(X, y)

        ## Main
        profiles = {}
        for class_i in np.sort(np.unique(y)):

            ## Subset for target label
            feature_matrix = X.loc[y == class_i, :]

            ## Scale high mut load samples down
            if self.count_ceiling is not None:
                feature_matrix = self._apply_count_ceiling(feature_matrix)

            ## Avg
            profile = self._calc_column_averages(feature_matrix)

            profiles[class_i] = profile

        profiles = pd.DataFrame.from_dict(profiles, orient="columns")
        profiles.index = X.columns

        if self.feature_prefix is not None:
            profiles.columns = self.feature_prefix + profiles.columns

        ## Normalize to sum to 1
        if self.normalize_profiles:
            profiles = profiles / np.sum(profiles, axis=0)

        self.profiles_ = profiles

        return self

    def transform(self, X: pd.DataFrame, y: Any = None) -> pd.DataFrame:
        """

        Parameters
        ----------
        X: pandas DataFrame of shape (n_samples, n_features)
            Where `n_samples` is the number of samples and `n_features` is the number of features

        y: None
            Not used. Argument only exists for compatibility

        Returns
        -------
        pandas DataFrame
            For each sample (row) similarities to each per cancer type profile (column)

        """
        check_is_fitted(self, "profiles_")
        X_trans = cos_sim(X, self.profiles_.transpose())
        return X_trans

    def fit_transform(self, X: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
        self.fit(X, y)
        return self.transform(X)

    def set_output(self, transform: Any = None) -> Self:
        return self

class NoiseProfileAdder(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        agg_func: Literal["sum","mean","median","iqm"] = "median",
        noise_counts: int | float = 100,
        count_ceiling: int | float = None,
        feature_prefix: str = None
    ):
        """
        Add background distribution of counts

        Inputs:
        - X, matrix of shape (n_samples, n_features) = feature matrix
        - y, array of shape (n_samples) = cancer type labels per sample

        At training time, compute:
        - profiles_per_class, matrix of shape (n_classes, n_features) = average mutation profile per cancer type in `y`
        - noise_profile, array of shape (n_features) = median of all per cancer type profiles = column medians of `X`

        At prediction time, add noise:
        - X_transformed = X + noise_counts * noise_profile

        Parameters
        ----------
        agg_func: str, 'sum', 'mean', 'iqm' (inter-quartile mean), or 'median'
            How to aggregate mutation profiles of the samples in each cancer type

        noise_counts: int
            Total number of noise mutations per sample to add

        count_ceiling: int
            When calculating the per cancer type profiles, the maximum total mutation count per sample. If a sample has
            more mutations than `count_ceiling`, the array will be scaled such that the sum of the array equals
            `count_ceiling`. This is so that hyper-mutator samples (i.e. those with extremely high number of mutations)
            do not dominate the cancer type profile.

        feature_prefix: str
            A prefix to append to the per cancer type profile matrix

        """
        self.agg_func = agg_func
        self.noise_counts = noise_counts
        self.count_ceiling = count_ceiling
        self.feature_prefix = feature_prefix

    def fit(self, X: pd.DataFrame, y: pd.Series) -> Self:

        transformer = ProfileSimilarityTransformer(
            agg_func=self.agg_func,
            count_ceiling=self.count_ceiling,
            normalize_profiles=True,
            feature_prefix=self.feature_prefix
        )
        profiles_per_class = transformer.fit(X, y).profiles_

        noise_profile = np.median(profiles_per_class, axis=1)
        self.noise_profile = noise_profile

        return self

    def transform(self, X: pd.DataFrame, y: Any = None):
        check_is_fitted(self, "noise_profile")

        X_trans = X
        X_trans = X_trans + self.noise_counts*self.noise_profile

        return X_trans

    def set_output(self, transform: str = None) -> Self:
        return self


class NonBestSimilarityScaler(BaseEstimator):
    def __init__(self, exponent: int | float = 1):
        """
        Scale rows of a matrix based on the max value of each row

        A major problem with cosine similarities is that they tend towards 1.0 the higher the number of dimensions
        (i.e. features). For gen_pos, the (sorted) cosine similarities look something like this

        >>> import numpy as np
        >>> cos_sims = np.array([0.95, 0.92, 0.87, 0.86, 0.75])

        We want the cosine similarities to not be so squished towards 1.0. In other words, we want to exaggerate the
        difference between consecutive cosine similarities. This improves the performance of the downstream logistic
        regression. To achieve this, we perform the following computations:

        >>> exponent = 5
        >>> max_value = np.max(cos_sims)
        >>> scale_factors = (cos_sims - max_value + 1)
        array([1., 0.97, 0.92, 0.91, 0.8])
        >>> scale_factors = scale_factors ** exponent
        array([1., 0.86, 0.67, 0.62, 0.32])
        >>> cos_sims_scaled = scale_factors * cos_sims
        array([0.95, 0.79 , 0.57, 0.54, 0.25])

        `exponent` is user defined, where are larger exponent scales the bottom values closer to zero

        The scaling results in the top values remaining high, while the bottom values are pushed closer to zero.

        Parameters
        ----------
        exponent : int or float
            where are larger exponent scales the bottom values closer to zero

        save_X_transformed : bool
            Store the transformation of `X` at `self.X_transformed`?

        """
        self.exponent = exponent

    def fit(self, X: Any, y: Any = None) -> Self:
        return self

    def transform(self, X: pd.DataFrame, y: Any = None):
        max_values = X.max(axis=1)
        scale_factors = X.add(-max_values + 1, axis=0)

        if self.exponent == 0:
            X_trans = X
        elif self.exponent == 1:
            X_trans = scale_factors * X
        else:
            X_trans = (scale_factors ** self.exponent) * X

        return X_trans

    def fit_transform(self, X: pd.DataFrame, y: Any = None, **fit_params) -> pd.DataFrame:
        return self.transform(X)

    def set_output(self, transform: Any = None) -> Self:
        return self


