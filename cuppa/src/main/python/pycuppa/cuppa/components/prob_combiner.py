import numpy as np
import pandas as pd
from typing import Self, Literal
from sklearn.base import BaseEstimator

from cuppa.constants import DEFAULT_FEATURE_PREFIX_SEPERATOR, META_CLF_NAMES
from cuppa.logger import LoggerMixin


class ProbCombiner(BaseEstimator, LoggerMixin):
    def __init__(
        self,
        combine_mode: Literal["multiply", "mean"] = "multiply",
        prob_floor: float = 0.01,
        sep: str = DEFAULT_FEATURE_PREFIX_SEPERATOR
    ):
        """

        Parameters
        ----------
        combine_mode : str
            "multiply": dna_combined and rna_combined probs below `prob_floor` are set to `prob_floor`, and then the
            probs are multiplied together. `prob_floor` serves to prevent multiplying zero probabilities, e.g. 0.99 * 0
            which results in the undesired combined prob of 0

            "mean": arithmetic mean of the dna_combined and rna_combined probs

        prob_floor : float
            See `combine_mode`

        sep : str, default = "__"
            How {prob_type}__{cancer_type} columns are separated

        """
        self.combine_mode = combine_mode
        self._check_combine_mode()

        self.prob_floor = prob_floor
        self.sep = sep

    def _check_combine_mode(self):
        if self.combine_mode not in ["multiply", "mean"]:
            self.logger.error("`combine_mode` must be 'multiply' or 'mean'")
            raise ValueError

    def _split_X_cols_by_prefix(self, X):
        affixes = X.columns.str.split(self.sep, n=1, expand=False)
        prefixes = affixes.map(lambda x: x[0])
        suffixes = affixes.map(lambda x: x[1])

        X_split = {}
        for prefix in prefixes.unique():
            selected_cols = prefixes==prefix

            X_subset = X.loc[:, selected_cols].copy()
            X_subset.columns = suffixes[selected_cols]

            X_split[prefix] = X_subset

        return X_split

    def transform(self, X, y=None):
        probs_split = self._split_X_cols_by_prefix(X)
        probs_dna = probs_split[META_CLF_NAMES.DNA_COMBINED]
        probs_rna = probs_split[META_CLF_NAMES.RNA_COMBINED]

        ## Select samples with RNA data
        na_rows = probs_rna.iloc[:, 0].isna()
        probs_rna = probs_rna.loc[~na_rows]
        probs_dna = probs_dna.loc[~na_rows, probs_rna.columns]
        ## Align columns as the RNA classifier may be missing cancer types present in the DNA classifier

        if self.prob_floor is not None:
            ## Account for cases where e.g. RNA prob is high, but DNA prob is 0
            probs_dna = probs_dna.copy()
            probs_rna = probs_rna.copy()

            probs_dna[probs_dna < self.prob_floor] = self.prob_floor
            probs_rna[probs_rna < self.prob_floor] = self.prob_floor

        if self.combine_mode == "multiply":
            probs_combined = probs_dna * probs_rna
        elif self.combine_mode == "mean":
            probs_combined = np.mean([probs_dna, probs_rna], axis=0)

        probs_combined = pd.DataFrame(
            probs_combined,
            index=probs_dna.index,
            columns=probs_dna.columns
        )

        ## Normalize 0 to 1
        row_sums = probs_combined.sum(axis=1, min_count=1)
        probs_combined = (probs_combined.T / row_sums).T

        return probs_combined

    def fit(self, X: pd.DataFrame, y: pd.Series = None) -> Self:
        return self

    def fit_transform(self, X: pd.DataFrame, y: pd.Series = None) -> Self:
        return self.transform(X)

    def set_output(self, transform: str = None) -> Self:
        return self
