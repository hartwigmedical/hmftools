from functools import cached_property
from typing import Self, Collection, Iterable, Optional
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from numpy.typing import NDArray
from cuppa.constants import DEFAULT_KEYWORDS_MALE_CLASSES, DEFAULT_KEYWORDS_FEMALE_CLASSES
from cuppa.logger import LoggerMixin


class FusionProbOverrider(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        sample_fusions: Optional[pd.DataFrame] = None,
        overrides: Optional[pd.DataFrame] = None,
        mask_base_value: float = 0.01,
        bypass: bool = False,

        _is_fitted: bool = False
    ):
        """
        When a cancer type specific fusion is found in a sample, diminish the probabilities of other cancer types

        The below example describes how probabilities are overriden.

        The RUNX1-RUNX1T1 fusion is characteristric of acute myeloid leukemia (AML). However, is often misclassified as
        myeloproliferative neoplasm (MPN) due to common genomic traits.

        Given a sample with a RUNX1-RUNX1T1 fusion, and the following prediction classes and probabilties:
        >>> classes = ["AML", "MPN", "Liver", "Prostate", ...]
        >>> probs = [0.26, 0.31, 0.02, 0.01, ...]

        Make a 'mask' array with:
            1.00 for target cancer type(s). Can be multiple cancer types
            0.01 for other cancer types. Value specified with the `mask_base_value` argument
        >>> mask = [1.00, 0.01, 0.01, 0.01, ...]

        Apply mask to the probabilties
        >>> probs = probs * mask ## [0.26, 0.03, 0.00, 0.00, ...]

        Make probs sum to 1
        >>> probs = probs / sum(probs)


        Parameters
        ----------
        sample_fusions: pandas DataFrame
            Fusions (columns) per sample (row). Each cell can be 1 (fusion present) or 0 (fusion absent).

        overrides: pandas DataFrame
            Specifies gene fusions and their target cancer types.
            Has the columns: feat_prefix, feat_basename, target_class

        mask_base_value: float
            See above description

        bypass: bool
            Do not override any probabilities?

        _is_fitted : bool
            This variable serves to silence the warning from _check_sample_fusions_empty() before cuppa classifier is
            trained
        """

        self.sample_fusions = sample_fusions

        self._is_fitted = _is_fitted
        self._check_sample_fusions_empty()

        self.overrides = overrides
        self._check_overrides_columns()
        self._check_overrides_empty()

        self.mask_base_value = mask_base_value
        self.bypass = bypass

    ## Check inputs ================================
    def _check_sample_fusions_empty(self) -> None:
        if self.sample_fusions is None and self._is_fitted:
            self.logger.debug("`sample_fusions` is None. Probs will not be overriden")

    def _check_overrides_empty(self) -> None:
        if self.overrides is None:
            self.logger.debug("`overrides` is None. Probs will not be overriden")

    def _check_overrides_columns(self) -> None:
        if self.overrides is not None:
            required_cols = pd.Series(["feat_prefix", "feat_basename", "target_class"])
            if not required_cols.isin(self.overrides).all():
                self.logger.error("`overrides` must have the columns: " + ", ".join(required_cols))
                raise KeyError

    ## Properties ================================
    @cached_property
    def _overrides_classes(self) -> pd.Series:
        return pd.Series(self.overrides["target_class"].sort_values().unique())

    def _check_missing_pred_classes(self, X):
        missing_pred_classes = self._overrides_classes[
            ~self._overrides_classes.isin(X.columns)
        ]

        if len(missing_pred_classes)>0:
            self.logger.warning("Some target classes were not found as pred classes: " + ", ".join(missing_pred_classes))

    @cached_property
    def _overrides_fusions(self) -> pd.Series:
        feat_names = self.overrides["feat_prefix"] + self.overrides["feat_basename"]
        return pd.Series(feat_names.sort_values().unique())

    def _get_existing_sample_fusions(self) -> pd.DataFrame:
        df = self.sample_fusions \
            .reindex(self._overrides_fusions, axis=1) \
            .stack().reset_index()
        df.columns = ["sample_id", "feat_name", "bool"]

        df = df[df["bool"] > 0]
        df = df.drop("bool", axis=1)

        return df

    ## Main ================================
    def transform(self, X: pd.DataFrame, y: Optional[pd.Series] = None) -> pd.DataFrame:
        # X = probs

        ## Init --------------------------------
        if self.overrides is None or self.sample_fusions is None or self.bypass:
            return X

        self._check_missing_pred_classes(X)

        existing_sample_fusions = self._get_existing_sample_fusions()
        if len(existing_sample_fusions)==0:
            return X

        ## Get which samples and pred classes to override
        overrides = pd.DataFrame(dict(
            feat_name = self.overrides["feat_prefix"] + self.overrides["feat_basename"],
            target_class = self.overrides["target_class"]
        ))

        override_info = existing_sample_fusions.merge(overrides, how="left", on="feat_name")

        ## Modify probabilities --------------------------------
        n_classes = X.shape[1]
        arr_template = np.empty(n_classes)
        arr_template.fill(self.mask_base_value)
        arr_template = pd.Series(arr_template, index=X.columns)

        ## Create output prob matrix. Convert multi-index to index
        X_trans = X.copy()
        if isinstance(X_trans.index, pd.MultiIndex):
            X_trans.index = X_trans.index.get_level_values("sample_id")

        ## Create mask array for each sample and multiply by the sample's feature array (i.e. row)
        sample_ids = override_info["sample_id"].unique()
        for sample_id in sample_ids:
            arr_mask = arr_template.copy()
            sample_target_classes = override_info.loc[override_info["sample_id"] == sample_id, "target_class"]
            arr_mask[sample_target_classes] = 1.0
            X_trans.loc[sample_id, :] *= arr_mask

        ## Make probs per sample sum to 1
        row_sums = X_trans.sum(axis=1, min_count=1)
        X_trans = (X_trans.T / row_sums).T

        ## Restore original (multi)-index
        X_trans.index = X.index

        return X_trans

    ## Required sklearn estimator methods ================================
    def fit(self, X: Optional[pd.DataFrame], y = None) -> Self:
        self._is_fitted = True
        return self

    def fit_transform(self, X: pd.DataFrame, y = None) -> Self:
        self._is_fitted = True
        return self.transform(X)

    def set_output(self, transform: None) -> Self:
        return self


class SexProbFilter(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        sample_sexes: Optional[pd.Series] = None,
        keywords_male_classes: Collection[str] = DEFAULT_KEYWORDS_MALE_CLASSES,
        keywords_female_classes: Collection[str] = DEFAULT_KEYWORDS_FEMALE_CLASSES,
        show_warnings: bool = False
    ):
        """
        Filter prediction probabilities by sample sex

        Suppose there is female sample, with a Prostate prediction with probability 0.1.

        As Prostate is a male cancer type, the probability of Prostate will be set to zero. This process is repeated for
        all male cancer type probabilities. The probabilities of this sample will then be re-normalized to sum to 1.

        Parameters
        ----------
        sample_sexes : pandas Series of type bool
            Sex of each sample, where True is male and False is female

        keywords_male_classes : array-like of type str
            Prediction classes containing these strings will be considered as male classes

        keywords_female_classes
            Prediction classes containing these strings will be considered as female classes

        """

        self.sample_sexes = sample_sexes
        self._check_sample_sexes_is_bool_series()

        self.keywords_male_classes = keywords_male_classes
        self.keywords_female_classes = keywords_female_classes

        self.show_warnings = show_warnings

    ## Checks ================================
    def _check_sample_sexes_empty(self):
        if self.sample_sexes is None and self.show_warnings:
            self.logger.debug("`sample_sexes` is None. Probs will not be overriden")

    def _check_sample_sexes_is_bool_series(self):
        #print(self.sample_sexes)
        if self.sample_sexes is None:
            return None

        if not isinstance(self.sample_sexes, pd.Series) or self.sample_sexes.dtype != "bool":
            self.logger.error("`sample_sexes` must be a pandas Series of type bool")
            raise TypeError

    ## Main ================================
    def _align_sample_sexes(self, X: pd.DataFrame) -> pd.Series:

        ## Align samples in `sample_sexes` and `X`
        sample_sexes_reindexed = self.sample_sexes.reindex(X.index) ## is np.nan when sample sex is missing for a given sample

        samples_missing_sex = sample_sexes_reindexed.isna().sum()
        if samples_missing_sex > 0 and self.show_warnings:
            self.logger.warning(
                "%s samples were missing in `sample_sexes.index` and did not have their probs overriden" %
                str(samples_missing_sex)
            )

        return sample_sexes_reindexed

    @staticmethod
    def _get_columns_matching_keywords(X: pd.DataFrame, keywords: Iterable) -> pd.Series:
        columns = X.columns
        keywords = pd.Series(keywords)

        if isinstance(columns, pd.MultiIndex):
            ## Look for matches at each multiindex level
            matches = [
                columns.get_level_values(level).str.contains("|".join(keywords), case=False)
                for level in range(columns.nlevels)
            ]
            matches = np.array(matches).any(axis=0)
        else:
            matches = columns.str.contains("|".join(keywords), case=False)

        return columns[matches]

    @staticmethod
    def _override_probs(X: pd.DataFrame, columns: Iterable, sex_mask: pd.Series | NDArray) -> pd.DataFrame:
        X_trans = X.copy() ## Prevent modifying the original copy

        ## If sex is unknown, leave probs unmodified; i.e. multiply by 1
        sex_mask = np.where(sex_mask.isna(), 1, sex_mask)

        for column in columns:
            X_trans[column] = X_trans[column] * sex_mask

        return X_trans

    @staticmethod
    def _normalize_rows_0_1(X: pd.DataFrame) -> pd.DataFrame:

        def divide_rows_by_rowsum(X):
            row_sums = X.sum(axis=1, min_count=1)
            X_trans = (X.T / row_sums).T
            return X_trans

        if isinstance(X.columns, pd.MultiIndex):
            ## Normalize by column group
            column_groups = X.columns.get_level_values(0).unique()
            X_trans = {
                group: divide_rows_by_rowsum(X[group])
                for group in column_groups
            }
            X_trans = pd.concat(X_trans, axis=1)

        else:
            X_trans = divide_rows_by_rowsum(X)

        return X_trans

    def transform(self, X: pd.DataFrame, y = None):

        sample_sexes_reindexed = self._align_sample_sexes(X)

        male_columns = self._get_columns_matching_keywords(X, self.keywords_male_classes)
        female_columns = self._get_columns_matching_keywords(X, self.keywords_female_classes)

        ## pd.Series of length n_samples. Is True when sample needs to have probs overriden
        male_mask = sample_sexes_reindexed
        female_mask = 1 - sample_sexes_reindexed

        X_trans = X.copy() ## Prevent modifying the original copy
        X_trans = self._override_probs(X_trans, male_columns, male_mask)
        X_trans = self._override_probs(X_trans, female_columns, female_mask)

        X_trans = self._normalize_rows_0_1(X_trans)

        return X_trans

    ## Required sklearn estimator methods ================================
    def fit(self, X: None, y: None) -> Self:
        return self

    def fit_transform(self, X: pd.DataFrame, y = None) -> Self:
        return self.transform(X)

    def set_output(self, transform: None) -> Self:
        return self




