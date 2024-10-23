from __future__ import annotations

import numpy as np
import pandas as pd
import sklearn.linear_model._logistic
from sklearn.utils.validation import check_is_fitted
from typing import Iterable
from sklearn.exceptions import ConvergenceWarning

import warnings
warnings.simplefilter("ignore", category=ConvergenceWarning)


class LogisticRegression(sklearn.linear_model._logistic.LogisticRegression):
    """
    Modified :class:`sklearn.linear_model.LogisticRegression`

    Changes:
    - :method:`predict_proba()` always returns a pandas DataFrame (default was 2d numpy array)
    - Added :method:`transform()` that mirrors :method:`predict_proba()`
    - Added :method:`feat_contrib()`.
    - :method:`fit()` now calculates the mean of each feature (required for :method:`feat_contrib`)

    See the below url for documentation of the original class:
    https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
    """

    def __init__(
        self,
        penalty: str = "l2",
        dual: bool = False,
        tol: float = 1e-4,
        C: float = 1.0,
        fit_intercept: bool = True,
        intercept_scaling: float = 1,
        class_weight: dict | str = None,
        random_state: int = None,
        solver: str = "lbfgs",
        max_iter: int = 100,
        multi_class: str = "auto",
        verbose: int = 0,
        warm_start: bool = False,
        n_jobs: int = None,
        l1_ratio: list[float] = None
    ):
        ## To inherit sklearn estimators, arguments need to explicitly be specified to the __init__ method
        ## https://stackoverflow.com/questions/40025406/inherit-from-scikit-learns-lassocv-model
        super(LogisticRegression, self).__init__(
            penalty = penalty,
            dual = dual,
            tol = tol,
            C = C,
            fit_intercept = fit_intercept,
            intercept_scaling = intercept_scaling,
            class_weight = class_weight,
            random_state = random_state,
            solver = solver,
            max_iter = max_iter,
            multi_class = multi_class,
            verbose = verbose,
            warm_start = warm_start,
            n_jobs = n_jobs,
            l1_ratio = l1_ratio
        )

    def fit(self, X: pd.DataFrame, y: pd.Series, sample_weight: Iterable[float] = None) -> "LogisticRegression":
        """
        Fit model

        Fits the coefficients and calculates the cohort-wide means for each feature

        Parameters
        ----------
        X: {array-like, sparse matrix} of shape (n_samples, n_features)
            Features (columns) for each sample (row)

        y: array-like of shape (n_samples,) default=None
            Labels for each sample

        sample_weight: array-like of shape (n_samples,) default=None
            Array of weights that are assigned to individual samples. If not provided, then each sample is given unit weight.

        Returns
        -------
        self
            Fitted estimator

        """
        super().fit(np.array(X), y, sample_weight)

        self.X_means = np.mean(X, axis=0)
        self.feature_names_in_ = np.array(X.columns)

        if len(self.classes_) == 2:
            ## Binary classification will only return a 1-row coefficient matrix which are for the positive class.
            ## Other methods expect no. of classes == no. rows in the coefficient matrix. Therefore, we need to make
            ## a 2 row matrix.
            coef_pos_class = self.coef_[0]
            coef_neg_class = coef_pos_class * -1
            self.coef_ = np.vstack([coef_neg_class, coef_pos_class])

        return self

    def predict_proba(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes.

        For a multi_class problem, if multi_class is set to be "multinomial"
        the softmax function is used to find the predicted probability of
        each class.
        Else use a one-vs-rest approach, i.e calculate the probability
        of each class assuming it to be positive using the logistic function.
        and normalize these values across all the classes.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Vector to be scored, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        pandas DataFrame of shape (n_samples, n_classes)
            Returns the probability of each sample for each class in the model

        """
        ## Handle dataframe with no rows
        if X.shape[0] == 0:
            return pd.DataFrame(index=None, columns=self.classes_, dtype="float64")

        if X.ndim == 1:
            index = X.index if isinstance(X, pd.DataFrame) else [0]

            ## Make 1-row matrix
            X = X.reshape(1, -1)
        else:
            index = X.index if isinstance(X, pd.DataFrame) else [i for i in range(X.shape[0])]

        probs = super().predict_proba(X)

        return pd.DataFrame(probs, index=index, columns=self.classes_)

    def transform(self, X: pd.DataFrame, y: pd.Series = None) -> pd.DataFrame:
        """

        Parameters
        ----------
        X: {array-like, sparse matrix} of shape (n_samples, n_features)
            Features (columns) for each sample (row)

        y: None
             Not used. Argument only exists for compatibility

        Returns
        -------
        pandas DataFrame of shape (n_samples, n_classes)
            Returns the probability of each sample for each class in the model

        """
        return self.predict_proba(X)

    def feat_imp(self, mode: str = "coef") -> pd.DataFrame:
        """
        Feature importance

        This method calculates SHAP values using the training data feature means using the following formula:
            coef[i] * X.mean(axis=0)[i] for the ith feature

            Where:
            - coef = coefficients of the logistic regression
            - X = feature matrix of the training data where rows are samples and columns are features. `X.mean(axis=0)`
            refers to the mean of each feature across all samples

        Reference: https://shap.readthedocs.io/en/latest/example_notebooks/overviews/An%20introduction%20to%20explainable%20AI%20with%20Shapley%20values.html

        Parameters
        ----------
        mode : str
            "shap": Shapley values
            "coef": Coefficients

        Returns
        -------
        pandas DataFrame of shape (n_classes, n_features)
            Where n_classes and n_features is the number of classes and features in the model respectively

        """

        check_is_fitted(self, "coef_")

        if mode not in ("shap", "coef"):
            raise ValueError("`mode` must be 'shap' or 'coef'")

        if mode=="shap":
            X_means = np.reshape(self.X_means, (1, -1))
            importance = self.coef_ * X_means
        else:
            importance = self.coef_

        importance = pd.DataFrame(importance, columns=self.feature_names_in_, index=self.classes_)

        return importance

    def feat_contrib(self, X: pd.DataFrame, y: pd.Series = None) -> pd.DataFrame:
        """

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

        y: None
             Not used. Argument only exists for compatibility

        Returns
        -------
        pandas multi-index DataFrame of shape (n_samples, n_classes) x n_features

        """

        check_is_fitted(self, "coef_")

        if not isinstance(X, pd.DataFrame):
            X = pd.DataFrame(X, columns=self.feature_names_in_)

        ## Get attributes from logistic regression
        n_classes = len(self.classes_)
        n_samples = X.shape[0]

        ## Feature contribution by prediction class
        contribs_by_class = []
        diff_vs_mean = X - self.X_means
        for class_index, class_name in enumerate(self.classes_):
            contribs_by_class.append(self.coef_[class_index] * diff_vs_mean)

        ## List of 2D arrays: n_samples (list items) x [n_pred_classes (rows), n_features (columns)]
        contribs_merged = np.array(pd.concat(contribs_by_class))
        del contribs_by_class

        sample_row_indexes = np.array([class_index * n_samples for class_index in range(0, n_classes)])
        contribs_by_sample = []
        for sample_i in range(n_samples):
            # contribs_sample_i = contribs_merged[sample_row_indexes]
            contribs_by_sample.append(contribs_merged[sample_row_indexes])
            sample_row_indexes = sample_row_indexes + 1
        del contribs_merged

        ## Multi indexed dataframe
        sample_names = np.repeat(X.index, n_classes)
        class_names = list(self.classes_) * n_samples
        indexes = pd.MultiIndex.from_arrays([sample_names, class_names])
        indexes.names = ["sample_id", "pred_class"]

        contribs = pd.DataFrame(
            np.vstack(contribs_by_sample),
            columns=self.feature_names_in_,
            index=indexes
        )

        ## Add priors
        priors = np.dot(self.coef_, self.X_means) + self.intercept_
        priors = np.array(list(priors) * n_samples)

        contribs.insert(0, "_prior", priors)

        ## Convert -0 to 0
        contribs += 0

        return contribs

    def set_output(self, transform: str = None) -> "LogisticRegression":
        return self
