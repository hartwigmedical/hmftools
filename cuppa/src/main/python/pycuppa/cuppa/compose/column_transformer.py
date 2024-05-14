from __future__ import annotations
from typing import Iterable

import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator
from sklearn.utils._estimator_html_repr import _VisualBlock
from sklearn.utils.metaestimators import _BaseComposition

from cuppa.constants import DEFAULT_FEATURE_PREFIX_SEPERATOR
from cuppa.logger import LoggerMixin


class ColumnTransformer(_BaseComposition, LoggerMixin):

    """
    Simplified ColumnTransformer based on the sklearn ColumnTransformer

    Applies transformers to columns of a pandas DataFrame.

    Parameters
    ----------
    transformers : list of tuples
        List of (name, transformer, columns) tuples specifying the
        transformer objects to be applied to subsets of the data.

        name : str
            Like in Pipeline and FeatureUnion, this allows the transformer and
            its parameters to be set using ``set_params`` and searched in grid
            search.
        transformer : {'drop', 'passthrough'} or estimator
            Estimator must support :term:`fit` and :term:`transform`.
            Special-cased strings 'drop' and 'passthrough' are accepted as
            well, to indicate to drop the columns or to pass them through
            untransformed, respectively.
        columns :  str, array-like of str, int, array-like of int, \
                array-like of bool, slice or callable
            Indexes the data on its second axis. Integers are interpreted as
            positional columns, while strings can reference DataFrame columns
            by name.  A scalar string or int should be used where
            ``transformer`` expects X to be a 1d array-like (vector),
            otherwise a 2d array will be passed to the transformer.
            A callable is passed the input data `X` and can return any of the
            above. To select multiple columns by name or dtype, you can use
            :obj:`make_column_selector`.

    n_jobs : int
        The number of jobs to use when calling `fit()`. One thread is used per transformer. If more threads are
        specified than the number of transformers, the number of transformers will be used as n_jobs instead

    verbose : bool, default=False
        Print progress messages?

    verbose_feature_names_out : bool, default=True
        If True, :meth:`ColumnTransformer.get_feature_names_out` will prefix
        all feature names with the name of the transformer that generated that
        feature.
        If False, :meth:`ColumnTransformer.get_feature_names_out` will not
        prefix any feature names and will error if feature names are not
        unique.

    prefix_sep : str, default="__".
        The separator to used when making output feature names in the format {column_group_name}__{feature_name},
        where column_group_name is `name` as provided in the `transformers` argument.

    """

    def __init__(
        self,
        transformers: list[tuple],
        n_jobs: int = 1,
        verbose: bool = False,
        verbose_feature_names_out: bool = True,
        prefix_sep = DEFAULT_FEATURE_PREFIX_SEPERATOR
    ):
        self.transformers = transformers
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.verbose_feature_names_out = verbose_feature_names_out
        self.prefix_sep = prefix_sep

    @property
    def transformer_names(self) -> list[str]:
        return [i[0] for i in self.transformers]

    def __getitem__(self, transformer_name):
        transformer_names = self.transformer_names

        if transformer_name not in transformer_names:
            self.logger.error("Invalid estimator name: %s. Valid names are: %s" % (
                transformer_name,
                ", ".join(transformer_names)
            ))
            raise ValueError

        for name, estimator, columns in self.transformers:
            if transformer_name == name:
                return estimator

    @staticmethod
    def _fit_one_transformer(
        X: pd.DataFrame,
        y: pd.Series,
        name: str,
        transformer: BaseEstimator,
        columns: Iterable[str],
        logger,
        verbose: bool
    ) -> tuple[str, BaseEstimator, Iterable[str]]:

        fitted_transformer = transformer.fit(X[columns], y)

        if verbose:
            logger.debug("Done fitting: " + name)

        return name, fitted_transformer, columns

    def _get_n_jobs_required(self, n_jobs: int) -> int:

        n_transformers = len(self.transformers)

        if n_jobs > n_transformers:
            self.logger.info("`n_jobs` (%i) is more than the number of transformers (%i). Using `n_jobs`=%i instead" % (
                n_jobs, n_transformers, n_transformers
            ))
            return n_transformers

        return n_jobs

    def fit(
        self,
        X: pd.DataFrame,
        y: pd.Series = None,
        n_jobs: int = None,
        verbose: bool = None
    ) -> "ColumnTransformer":
        """

        Parameters
        ----------
        X: pandas DataFrame
            Features (columns) per sample (row)

        y: pandas Series
            Sample labels

        n_jobs: int
            Number of threads to use. One thread is used per transformer. If more threads are specified than the
            number of transformers, the number of transformers will be used as n_jobs instead

        verbose: bool
            Print progress messages?

        Returns
        -------
        self

        """
        if verbose is None:
            verbose = self.verbose

        if verbose:
            self.logger.info("Fitting transformers: " + ", ".join([name for name,_,_ in self.transformers]))

        self.feature_names_in_ = X.columns

        if n_jobs is None:
            n_jobs = self.n_jobs
        n_jobs = self._get_n_jobs_required(n_jobs)

        if n_jobs == 1:
            for name, transformer, columns in self.transformers:
                self._fit_one_transformer(
                    X=X, y=y,
                    name=name, transformer=transformer, columns=columns,
                    logger=self.logger, verbose=verbose
                )
            return self

        parallel = Parallel(
            backend="threading",
            n_jobs=n_jobs,
            verbose=self.verbose
        )

        fitted_transformers = parallel(
            delayed(self._fit_one_transformer)(
                X=X, y=y,
                name=name, transformer=transformer, columns=columns,
                logger=self.logger,
                verbose=verbose
            )
            for name, transformer, columns in self.transformers
        )

        self.transformers = fitted_transformers

    def transform(
        self,
        X: pd.DataFrame,
        y: pd.Series = None,
        return_Xs: bool = False,
        verbose_feature_names_out: bool = None,
        verbose: bool = None
    ) -> pd.DataFrame | dict[str, pd.DataFrame]:

        """

        Parameters
        ----------
        X: pandas DataFrame
            Features (columns) per sample (row)

        y: pandas Series
            Sample labels

        return_Xs: bool
            Return a dataframe for each column group transformation (as a dict)

        verbose_feature_names_out : bool, default=True
            If True, :meth:`ColumnTransformer.get_feature_names_out` will prefix
            all feature names with the name of the transformer that generated that
            feature.
            If False, :meth:`ColumnTransformer.get_feature_names_out` will not
            prefix any feature names and will error if feature names are not
        unique.

        verbose: bool
            Print progress messages?

        Returns
        -------
        Pandas dataframe or dict of pandas dataframes

        """

        if verbose is None:
            verbose = self.verbose

        if verbose_feature_names_out is None:
            verbose_feature_names_out = self.verbose_feature_names_out

        if verbose:
            self.logger.info("Processing transformers: " + ", ".join([name for name, _, _ in self.transformers]))

        Xs = {}

        for name, transformer, columns in self.transformers:
            if verbose:
                self.logger.debug("Transforming: " + name)

            X_trans = transformer.transform(X[columns])

            if verbose_feature_names_out:
                X_trans.columns = name + self.prefix_sep + X_trans.columns

            Xs[name] = X_trans

        if return_Xs:
            return Xs

        return pd.concat(Xs.values(), axis="columns")


    def fit_transform(self, X: pd.DataFrame, y: pd.Series = None) -> pd.DataFrame:
        self.fit(X, y)
        return self.transform(X)

    def feat_contrib(self, X: pd.DataFrame, show_warnings: bool = False) -> pd.DataFrame | None:
        """
        Feature contributions from transformers

        For each transformer, the feat_contrib() method is called if it exists for that transformer.

        Parameters
        ----------
        X: pandas DataFrame
            Features (columns) per sample (row)

        show_warnings:
            Show warning when no transformers have the feat_contrib() method?

        Returns
        -------

        """
        feat_contribs = {}
        for name, estimator, columns in self.transformers:
            if hasattr(estimator, "feat_contrib"):
                feat_contribs_i = estimator.feat_contrib(X[columns])

                ## Pipeline objects without any estimators with the feat_contrib() method return None
                ## This is so that nested Pipeline.feat_contrib() calls do not crash outer Pipeline.feat_contrib() calls
                if feat_contribs_i is None:
                    continue

                feat_contribs[name] = feat_contribs_i

        if len(feat_contribs)==0:
            if show_warnings:
                self.logger.warn("No estimators in `self.transformers` were found with the method `feat_contrib()`. Returning `None`")
            ## Also return None so that outer Pipeline.feat_contrib() calls do not crash
            return None

        feat_contribs = pd.concat(feat_contribs, axis=1)
        return feat_contribs

    def _sk_visual_block_(self) -> _VisualBlock:
        names, transformers, name_details = zip(*self.transformers)
        return _VisualBlock("parallel", transformers, names=names, name_details=name_details)

    def set_output(self, transform = None) -> "ColumnTransformer":
        return self


