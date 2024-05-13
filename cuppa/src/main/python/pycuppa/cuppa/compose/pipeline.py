from __future__ import annotations

import threading

import pandas as pd
import logging
import os
import shutil
from functools import cached_property
from typing import Iterable, Optional

import sklearn.pipeline
from sklearn import clone
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection._split import _BaseKFold
from sklearn.utils.parallel import delayed, Parallel

from cuppa.classifier.cuppa_prediction import CuppaPrediction
from cuppa.logger import LoggerMixin
import joblib

logger = logging.getLogger(__name__)

## Pipeline ================================
class Pipeline(sklearn.pipeline.Pipeline, LoggerMixin):
    """
    Modified :class:`sklearn.pipeline.Pipeline`

    Additions:
    - Option to cache pipeline steps so that they can be reloaded. Useful when certain steps are compute heavy
    - Option to fit and/or transform both `X` and `y`. In the original class, only `X` could be fitted/transformed
    - :method:`transform()` and :method:`fit_transform()` can transform until a specified step and/or keep intermediate
    transformed features

    See the below url for documentation of the original class:
    https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html


    Added parameters
    ----------------
    fit_cache_dir : str
        Path to a directory for caching during fitting. The directory will be created if it doesn't exist (as long as
        parent directories exist)
    """

    def __init__(
        self,
        steps: list[tuple],
        *,
        fit_cache_dir: Optional[str] = None,
        verbose: bool = False
    ):
        super().__init__(steps=steps, verbose=verbose)

        ## Fit cache
        self.fit_cache_dir = fit_cache_dir
        if self.fit_cache_dir is not None:
            os.makedirs(self.fit_cache_dir, exist_ok=True)

        self.set_output(transform="pandas")

    @staticmethod
    def fit_cached(
        estimator,
        X: pd.DataFrame,
        y: pd.Series,
        cache_path: Optional[str] = None,
        step_name: str = "estimator",
        logger = logger,
        verbose: bool = True,
        **fit_params
    ) -> "Pipeline":
        cache_exists = cache_path is not None and os.path.exists(cache_path)

        ## Load cache
        if cache_exists:
            estimator = joblib.load(cache_path)
            if verbose:
                logger.info("Loaded cache: " + step_name)
            return estimator

        ## Fit and write cache
        estimator = clone(estimator)  ## Clone to ensure estimator can be pickled
        estimator.fit(X, y, **fit_params)

        if cache_path is not None:
            joblib.dump(estimator, cache_path, compress=3)

        if verbose:
            logger.debug("Done fitting: " + step_name)

        return estimator

    def fit(self, X: pd.DataFrame, y: pd.Series, verbose: bool = None) -> "Pipeline":

        if verbose is None:
            verbose = self.verbose

        if verbose:
            self.logger.info("Fitting transformers: " + ", ".join([step_name for step_name,_ in self.steps]))

        self.steps = list(self.steps)
        self._validate_steps()

        ## Get **fit_params for each step
        fit_params_steps = {name: {} for name, step in self.steps if step is not None}

        ## Get the character width of the total number of steps
        n_steps = len(self.steps)
        n_steps_nchar = len(str(n_steps))

        ## Fit
        X_trans = X
        y_trans = y
        for step_idx, (step_name, estimator) in enumerate(self.steps):
            if estimator is None or estimator == "passthrough":
                continue

            ## Set up / load cache
            cache_path = None
            if self.fit_cache_dir is not None:
                if not os.path.exists(self.fit_cache_dir):
                    os.mkdir(self.fit_cache_dir)

                step_idx_string = str(step_idx + 1).rjust(n_steps_nchar, "0")
                cache_path = os.path.join(self.fit_cache_dir, step_idx_string + "-" + step_name + ".pickle.gz")

            ## Fit estimator
            fitted_estimator = self.fit_cached(
                estimator=estimator,
                X=X_trans,
                y=y_trans,
                cache_path=cache_path,

                step_name=step_name,
                verbose=verbose,
                logger=self.logger,

                **fit_params_steps[step_name],
            )

            ## Overwrite X (and y) with transformed version
            if hasattr(fitted_estimator, "fit_resample"):
                X_trans = fitted_estimator.fit_resample(X,y)
            else:
                X_trans = fitted_estimator.transform(X_trans, y_trans)

            ## Unpack if X is a tuple, and overwrite y
            if isinstance(X_trans, tuple):
                X_trans, y_trans = X_trans

            self.steps[step_idx] = (step_name, fitted_estimator)

        return self

    def clear_fit_cache(self) -> None:
        if self.fit_cache_dir is None:
            self.logger.warning("`self.fit_cache_dir` was not provided")
            return None

        os.rmdir(self.fit_cache_dir)

    @staticmethod
    def check_step_names(
        steps: Iterable[str],
        valid_steps: Iterable[str],
        logger = logger
    ) -> None:
        # steps = ["sub_clfs", "post_process", "fashkjdfa"]
        # valid_steps = self.named_steps.keys()

        steps = pd.Series(steps)
        is_valid = steps.isin(valid_steps)

        if any(~is_valid):
            invalid_steps = steps[~is_valid]

            error_msg = "The following keys are invalid: %s. Valid keys are: %s" % (
                ', '.join(f"'%s'" % w for w in invalid_steps),
                ', '.join(f"'%s'" % w for w in valid_steps)
            )
            logger.error(error_msg)
            raise KeyError

    def transform(
        self,
        X: pd.DataFrame,
        y: Optional[pd.Series] = None,
        return_y: bool = False,
        until_step: Optional[str] = None,
        keep_steps: Optional[str | list[str]] = None,
        verbose: Optional[bool] = None
    ) -> pd.DataFrame | tuple[pd.DataFrame, pd.Series]:
        """
        Parameters
        ----------
        X: pandas DataFrame
            Features (columns) per sample (row)

        y: pandas Series
            Sample labels

        return_y: bool
            Return transformed `y`? If True, a tuple with transformed X and y is returned

        until_step: str
            Step name to transform until

        keep_steps: str or list[str]
            Which steps to keep? If not None, will return a dict of data frames of transformed features at each step

        verbose: bool
            Show progress messages?

        """
        if verbose is None:
            verbose=self.verbose

        if until_step is not None:
            self.check_step_names(
                steps = until_step,
                valid_steps = self.named_steps.keys(),
                logger = self.logger
            )

        Xs = {}
        if keep_steps is not None:
            self.check_step_names(
                steps = keep_steps,
                valid_steps = self.named_steps.keys(),
                logger = self.logger
            )

        ## Main
        for index, name, estimator in self._iter():

            if verbose:
                self.logger.info("Transforming: " + name)

            X = estimator.transform(X, y)

            if verbose:
                self.logger.info("Transforming: " + name + " (done)")

            if keep_steps is not None and name in keep_steps:
                Xs[name] = X

            ## Unpack X,y if y is transformed, and overwrite y
            if isinstance(X, tuple):
                X, y = X

            ## Stop transform prematurely
            if until_step is not None and name==until_step:
                break

        ## Output
        if keep_steps is not None:
            return Xs

        if return_y:
            return X, y
        return X

    def fit_transform(
        self,
        X: pd.DataFrame,
        y: Optional[pd.Series] = None,
        return_y: bool = False,
        until_step: Optional[str | list[str]] = None,
        **fit_params
    ) -> pd.DataFrame:
        if until_step is not None:
            self.check_step_names(
                steps = until_step,
                valid_steps = self.named_steps.keys(),
                logger = self.logger
            )

        self.fit(X, y, **fit_params)
        return self.transform(X, y, return_y=return_y, until_step=until_step)

    def predict_proba(self, X: pd.DataFrame, **predict_proba_params) -> pd.DataFrame:
        X_trans = X

        for index, name, estimator in self._iter(with_final=False):
            X_trans = estimator.transform(X_trans)

        return self.steps[-1][1].predict_proba(X_trans, **predict_proba_params)

    def feat_contrib(self, X: pd.DataFrame, y: Optional[pd.Series] = None) -> pd.DataFrame | None:

        ## Search for estimators with `feat_contrib()` method
        has_feat_contrib_method = False
        for name, estimator in reversed(self.steps):
            if hasattr(estimator, "feat_contrib"):
                has_feat_contrib_method = True
                break

        if not has_feat_contrib_method:
            # if show_warnings:
            #     _LOGGER.log(
            #         "No step in `self.steps` was found with the method `feat_contrib()`. Returning `None`",
            #         level="warn"
            #     )

            ## Return None so that nested Pipeline.feat_contrib() calls do not crash outer Pipeline.feat_contrib() calls
            return None

        ## Get step name before the estimator with the `feat_contrib()` method and transform X until that step
        step_names = pd.Series(self.named_steps.keys())
        previous_step_index = step_names.index[step_names==name] - 1 ## This equals -1 when the Pipeline has exactly 1 estimator with the `feat_contrib()` method
        previous_step_index = previous_step_index[0]
        X_trans = self.transform(X, until_step=step_names[previous_step_index]) if previous_step_index != -1 else X

        ##
        contribs = estimator.feat_contrib(X_trans)
        return contribs

    def _get_estimators(self) -> list[str]:
        return [estimator for name, estimator in self.steps]

    def _get_all_estimators(self) -> list[str]:
        estimator_list = []

        estimator_list += self._get_estimators()

        for estimator in estimator_list:
            if hasattr(estimator, "_get_estimators"):
                estimator_list += estimator._get_estimators()

        return estimator_list


## CV ================================
class PipelineCrossValidator(LoggerMixin):
    def __init__(
        self,
        pipeline: Pipeline,
        X: pd.DataFrame,
        y: pd.Series,
        out_dir: Optional[str] = None,
        cv: _BaseKFold = StratifiedKFold(n_splits=10),
        y_split: Optional[pd.Series] = None,
        n_jobs: int = -1,
        pre_dispatch: str ="2*n_jobs",
        verbose: bool = True
    ):
        self.pipeline = pipeline
        self.X = X
        self.y = y
        self.out_dir = out_dir

        self.cv = cv
        self.y_split = y if y_split is None else y_split

        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self.verbose = verbose

    @cached_property
    def _cv_splits(self) -> list[int]:
        return list(self.cv.split(self.X, self.y_split))

    @cached_property
    def _n_cv_folds(self) -> int:
        return self.cv.get_n_splits()

    def get_one_test_set(self, cv_index: int) -> tuple[pd.DataFrame, pd.Series]:
        test_indexes = self._cv_splits[cv_index][1]
        X = self.X.iloc[test_indexes]
        y = self.y.iloc[test_indexes]

        return X,y

    def get_one_train_set(self, cv_index: int) -> tuple[pd.DataFrame, pd.Series]:
        train_indexes = self._cv_splits[cv_index][0]
        X = self.X.iloc[train_indexes]
        y = self.y.iloc[train_indexes]

        return X, y

    def _fit_one_pipeline(self, cv_index, **fit_params) -> "Pipeline":

        threading.current_thread().name = "CvFold-" + str(cv_index+1)

        X_train, y_train = self.get_one_train_set(cv_index)
        pipeline = self.pipelines[cv_index]
        return pipeline.fit(X_train, y_train, **fit_params)

    def fit(self, cache_training: bool = True, overwrite: bool = False, **fit_params) -> None:

        ## Output dir --------------------------------
        if self.verbose:
            self.logger.info("Training data: %s samples, %s features, %s classes" % (
                str(self.X.shape[0]),
                str(self.X.shape[1]),
                len(self.y.unique())
            ))

        if self.out_dir is not None:
            self.logger.info("Outputting to: " + self.out_dir)

            if overwrite and os.path.exists(self.out_dir):
                shutil.rmtree(self.out_dir, ignore_errors=True)

            if not os.path.exists(self.out_dir):
                os.makedirs(self.out_dir, exist_ok=True)

        ## Clone pipelines --------------------------------
        n_cv_folds = self.cv.get_n_splits()

        self.pipelines = []
        for cv_index in range(0, n_cv_folds):
            pipeline = clone(self.pipeline)
            self.pipelines.append(pipeline)

        ## Set up variables for saving CV models --------------------------------
        if cache_training:
            ## Parent dir
            model_cache_dir = os.path.join(self.out_dir, "models_cache/")
            if not os.path.exists(model_cache_dir):
                os.mkdir(model_cache_dir)

            ##
            max_cv_index_nchar = len(str(n_cv_folds + 1))
            for cv_index, pipeline in enumerate(self.pipelines):
                # cv_index=0
                ## Clone pipeline and assign cache dir
                cv_fold_name = str(cv_index+1).rjust(max_cv_index_nchar, "0")
                fit_cache_dir = os.path.join(model_cache_dir, cv_fold_name)
                pipeline.fit_cache_dir = fit_cache_dir

                if not os.path.exists(fit_cache_dir):
                    os.mkdir(fit_cache_dir)

        ## Main --------------------------------
        parallel = Parallel(
            backend="threading",
            n_jobs=self.n_jobs,
            pre_dispatch=self.pre_dispatch,
            verbose=False
        )

        pipelines = parallel(
            delayed(self._fit_one_pipeline)(cv_index=cv_index, **fit_params)
            for cv_index in range(0, self._n_cv_folds)
        )

        self.cv_models = pipelines

    def _check_is_cross_validated(self) -> None:
        if not hasattr(self, "cv_models"):
            self.logger.error("Attribute `cv_models` does not exist, likely because `cross_validate()` has not yet been called")
            raise AttributeError

    def _apply_on_one_test_set(self, method_name: str, cv_index: int, *args, **kwargs) -> pd.DataFrame:

        X_test, y_test = self.get_one_test_set(cv_index)
        cv_model = self.cv_models[cv_index]

        if self.verbose:
            test_set_num = cv_index + 1
            self.logger.info(f"Applying method `{method_name}` to test set {test_set_num}")

        method = getattr(cv_model, method_name)
        return method(X=X_test, y=y_test, *args, **kwargs)

    def apply_on_test_sets(
        self,
        method_name: str,
        concat_output: bool = True,
        n_jobs: int = -1,
        pre_dispatch: str = "2*n_jobs",
        verbose: bool = True,
        *args, **kwargs
    ) -> pd.DataFrame | CuppaPrediction:
        self._check_is_cross_validated()

        if not hasattr(self.pipeline, method_name):
            self.logger.error("`%s` not found as a method in `self.pipeline`" % method_name)
            raise AttributeError

        ## Iterate over test sets --------------------------------
        parallel = Parallel(n_jobs=n_jobs, pre_dispatch=pre_dispatch, verbose=verbose)
        X_y_trans = parallel(
            delayed(self._apply_on_one_test_set)(method_name=method_name, cv_index=cv_index, *args, **kwargs)
            for cv_index in range(0, self._n_cv_folds)
        )

        ## Concatenate output --------------------------------
        def _concat_X(X_list) -> pd.DataFrame | CuppaPrediction:
            #X_list=X_y_trans
            # X_list[0].__class__ == CuppaPrediction
            # type(X_list[0]) is cuppa.cuppa_prediction.CuppaPrediction
            # type(cuppa.cuppa_prediction.CuppaPrediction)

            if isinstance(X_list[0], CuppaPrediction):
                concat = CuppaPrediction.concat
            else:
                concat = pd.concat

            return concat(X_list)

        def _concat_y(y_list: list[pd.Series]) -> pd.Series:
            return pd.concat(y_list)

        if concat_output:
            if isinstance(X_y_trans[0], tuple):
                X_trans = _concat_X([i[0] for i in X_y_trans])
                y_trans = _concat_y([i[1] for i in X_y_trans])

                output = (X_trans, y_trans)
            else:
                output = _concat_X(X_y_trans)

        else:
            output = X_y_trans

        return output

    def transform(
        self,
        return_y: bool = False,
        until_step: Optional[str] = None,
        keep_steps: Optional[str | list[str]] = None
    ) -> pd.DataFrame:
        X_trans = self.apply_on_test_sets(
            "transform",
            return_y=return_y, until_step=until_step, keep_steps=keep_steps,
            concat_output=True, n_jobs=self.n_jobs, pre_dispatch=self.pre_dispatch, verbose=self.verbose
        )

        return X_trans

    def predict_proba(self) -> pd.DataFrame:

        probs = self.apply_on_test_sets(
            "predict_proba",
            concat_output=True, n_jobs=self.n_jobs, pre_dispatch=self.pre_dispatch, verbose=self.verbose
        )

        return probs

    def feat_imp(self, mode: str = "coef", return_means: bool = False) -> pd.DataFrame:

        ## Get feature importances across CV folds
        importances = [cv_model.feat_imp(mode=mode) for cv_model in self.cv_models]

        ## Convert to multi index dataframe
        importances = dict(enumerate(importances, start=1))
        importances = pd.concat(importances)
        importances.index.names = ("cv_fold", "class")

        ## Mean across CV folds
        if return_means:
            importances = importances.fillna(0).groupby("class").mean()

        return importances

    def get_cal_curves(self) -> pd.DataFrame:
        cal_curves = []
        for cv_index, model in enumerate(self.cv_models):
            cal_curves_i = model.get_cal_curves()
            cal_curves_i.insert(0, "group", cv_index)
            cal_curves.append(cal_curves_i)

        cal_curves = pd.concat(cal_curves)
        return cal_curves


