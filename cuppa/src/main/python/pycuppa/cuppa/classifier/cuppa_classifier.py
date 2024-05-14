from __future__ import annotations

import os
from typing import TYPE_CHECKING, Optional, Literal

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import make_column_selector

import cuppa.compose.pipeline
from cuppa.constants import SUB_CLF_NAMES, META_CLF_NAMES, LAYER_NAMES, CLF_GROUPS, NA_FILL_VALUE, \
    SIG_QUANTILE_TRANSFORMER_NAME, CLF_NAMES
from cuppa.classifier.classifiers import ClassifierLayers, SubClassifiers, MetaClassifiers
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredictionBuilder
from cuppa.classifier.feature_importance import FeatureImportance
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.components.mutational_signatures import SigCohortQuantileTransformer
from cuppa.components.prob_overriders import FusionProbOverrider, SexProbFilter
from cuppa.compose.column_transformer import DEFAULT_FEATURE_PREFIX_SEPERATOR
from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH, SEX_FEATURE_NAME
from cuppa.classifier.cuppa_classifier_utils import MissingFeaturesHandler, BypassedClassifierBuilder

if TYPE_CHECKING:
    from cuppa.performance.performance_stats import PerformanceStats


class CuppaClassifier(cuppa.compose.pipeline.Pipeline):
    """
    CUPPA classifier

    The CUPPA classifier is composed of 3 main layers of multinomial logistic regressions.

    Layer 1
        This layer contains 'sub-classifiers' that output a preliminary probability for each feature type:
            DNA
            - gen_pos: genomic position. SNV counts in each of the ~6000 x 500kb bins across the genome
            - snv96: counts of SNVs in each of the 96 trinucleotide contexts
            - event: various mutation events, including driver mutations, gene fusions, viral sequence insertions,
            structural variant (SV) types, tumor mutational burden (TMB), whole genome duplication, and sample sex

            RNA
            - gene_exp: gene expression. log(TPM+1) for each gene (TPM: transcripts per million)
            - alt_sj: alt splice junction. log(number of reads + 1)

        Each sub-classifier has the following general structure:
        - Feature transformations (e.g. scaling, log transform, etc.)
        - Logistic regression with lasso (aka L1) regularization

    Layer 2
        2 'meta-classifiers':
        - a DNA classifier that weighs the importance of the 3 DNA sub-classifiers
        - an RNA classifier that weighs the importance of the 2 RNA sub-classifiers

        Each meta-classifier has the following general structure:
        - Logistic regression with ridge (aka L2) regularization
        - Probability calibration
        - Sex filter

        The DNA classifier also has a fusion overrider which boosts the probability of e.g. Acute myeloid leukemia
        if a RUNX1-RUNXT1 fusion is present in the respective sample

    Layer 3
        A `ProbCombiner` that multiplies the probabilities of the DNA and RNA classifiers to yield a combined
        probability

    Parameters
    ----------
    fusion_overrides_path : str or None
        Path to the fusion overrides file containing the columns: feat_prefix, feat_basename, target_class
        If no path is specified, fusion override will be bypassed

    fit_cache_dir : str
        Path to a directory for caching. The directory will be created if it doesn't exist (as long as parent
        directories exist)

    n_jobs : int
        Number of threads to use when fitting

    verbose : bool
        Print progress messages when calling fit?

    Returns
    -------
    self
        The CUPPA classifier object, which is a child of the `sklearn.pipeline.Pipeline` class

    """

    def __init__(
        self,
        fusion_overrides_path: Optional[str] = DEFAULT_FUSION_OVERRIDES_PATH,
        fit_cache_dir: Optional[str] = None,
        n_jobs: int = 1,
        verbose: bool = True
    ):
        self.sig_quantile_transformer = None
        self.fusion_overrides_path = fusion_overrides_path
        self.n_jobs = n_jobs
        self.verbose = verbose

        self.cv_performance: PerformanceStats | None = None

        steps = [
            (LAYER_NAMES.SUB_CLF,
             ClassifierLayers.SubClassifierLayer(n_jobs=n_jobs, verbose=self.verbose)),

            (LAYER_NAMES.META_CLF,
             ClassifierLayers.MetaClassifierLayer(fusion_overrides_path=self.fusion_overrides_path, n_jobs=n_jobs, verbose=self.verbose)),

            (LAYER_NAMES.COMBINED, ClassifierLayers.CombinedLayer())
        ]

        super().__init__(steps=steps, fit_cache_dir=fit_cache_dir, verbose=verbose)

    def __repr__(self):
        return f"{self.__class__.__name__}(steps={self.steps.__repr__()})"

    ## I/O ================================
    def to_file(self, path: str) -> None:
        if self.verbose:
            self.logger.info("Exporting classifier to: " + path)
        joblib.dump(self, filename=path, compress=9)

    @staticmethod
    def _check_is_pickle_file(path: str) -> None:
        if not path.endswith(".pickle") and not path.endswith(".pickle.gz"):
            raise ValueError("Path must end with .pickle or .pickle.gz")

    @classmethod
    def from_file(cls, path: str) -> "CuppaClassifier":
        cls.get_class_logger(cls).info("Loading classifier from: " + str(path))
        cls._check_is_pickle_file(path)
        return joblib.load(path)

    ## Classifier accessors ================================
    @property
    def gen_pos_clf(self) -> SubClassifiers.GenPosClassifier:
        return self[LAYER_NAMES.SUB_CLF][SUB_CLF_NAMES.GEN_POS]

    @property
    def snv96_clf(self) -> SubClassifiers.Snv96Classifier:
        return self[LAYER_NAMES.SUB_CLF][SUB_CLF_NAMES.SNV96]

    @property
    def event_clf(self) -> SubClassifiers.EventClassifier:
        return self[LAYER_NAMES.SUB_CLF][SUB_CLF_NAMES.EVENT]

    @property
    def gene_exp_clf(self) -> SubClassifiers.GeneExpClassifier:
        return self[LAYER_NAMES.SUB_CLF][SUB_CLF_NAMES.GENE_EXP]

    @property
    def alt_sj_clf(self) -> SubClassifiers.AltSjClassifier:
        return self[LAYER_NAMES.SUB_CLF][SUB_CLF_NAMES.ALT_SJ]

    @property
    def sub_clf_layer(self) -> ClassifierLayers.SubClassifierLayer:
        return self[LAYER_NAMES.SUB_CLF]

    @property
    def dna_combined_clf(self) -> MetaClassifiers.DnaCombinedClassifier:
        return self[LAYER_NAMES.META_CLF][META_CLF_NAMES.DNA_COMBINED]

    @property
    def rna_combined_clf(self) -> MetaClassifiers.RnaCombinedClassifier:
        return self[LAYER_NAMES.META_CLF][META_CLF_NAMES.RNA_COMBINED]

    @property
    def meta_clf_layer(self) -> ClassifierLayers.MetaClassifierLayer:
        return self[LAYER_NAMES.META_CLF]

    @property
    def combined_layer(self) -> ClassifierLayers.CombinedLayer:
        return self[LAYER_NAMES.COMBINED]

    ## Overrides ================================
    ## Sex filter --------------------------------
    @property
    def sex_filters(self) -> list[SexProbFilter]:
        return [
            self.dna_combined_clf["sex_filter"],
            self.rna_combined_clf["sex_filter"]
        ]

    def set_sample_sexes(self, X: pd.DataFrame) -> None:
        sample_sexes = X[SEX_FEATURE_NAME].astype(bool)

        for sex_filter in self.sex_filters:
            sex_filter.sample_sexes = sample_sexes

    def reset_sample_sexes(self) -> None:
        for sex_filter in self.sex_filters:
            sex_filter.sample_sexes = None

    ## Fusions --------------------------------
    @property
    def fusion_prob_overrider(self) -> FusionProbOverrider:
        return self.dna_combined_clf["fusion_overrider"]

    def set_sample_fusions(self, X: pd.DataFrame) -> None:
        sample_fusions = X[make_column_selector("^event[.]fusion[.]")]
        self.fusion_prob_overrider.sample_fusions = sample_fusions

    def reset_sample_fusions(self) -> None:
        self.fusion_prob_overrider.sample_fusions = None

    ## Calibrators --------------------------------
    @property
    def prob_calibrators(self) -> dict[str, RollingAvgCalibration]:
        return dict(
            dna_combined = self.dna_combined_clf["calibrator"],
            rna_combined = self.rna_combined_clf["calibrator"]
        )

    def get_cal_curves(self) -> pd.DataFrame:
        calibrators = self.prob_calibrators

        cal_curves = []
        for clf_name, calibrator in calibrators.items():
            cal_curves_i = calibrator.get_cal_curves()
            cal_curves_i.insert(0,"clf_name",clf_name)
            cal_curves.append(cal_curves_i)

        cal_curves = pd.concat(cal_curves)

        return cal_curves

    ## Training ================================
    def classes_(self) -> None:
        ## `classes_` is a read-only @property in sklearn.pipeline.Pipeline. Undo this so that it can be overriden
        pass

    def fit_sig_quantile_transformer(self, X: pd.DataFrame, y: pd.Series) -> "SigCohortQuantileTransformer":
        ## Mutational signature quantiles per class
        ## These are not used for prediction and therefore are generated outside the Pipeline.fit() call
        transformer = SigCohortQuantileTransformer(clip_upper=False, clip_lower=False)

        transformer = cuppa.compose.pipeline.Pipeline.fit_cached(
            estimator=transformer,
            X=X[make_column_selector("^sig")],
            y=y,
            cache_path=None if self.fit_cache_dir is None else os.path.join(self.fit_cache_dir, "0-sig_quantile_transformer.pickle.gz"),
            verbose=self.verbose,
            step_name="sig_quantile_transformer",
            logger=self.logger
        )

        return transformer

    def fit(self, X: pd.DataFrame, y: pd.Series) -> "CuppaClassifier":

        if self.verbose:
            self.logger.info("Training CuppaClassifier on: %s samples, %s features, %s classes" % (
                str(X.shape[0]),
                str(X.shape[1]),
                len(y.unique())
            ))

            self.logger.info(
                "Fitting transformers: sig_quantile_transformer, " +
                ", ".join([step_name for step_name, _ in self.steps])
            )

        self.classes_ = np.unique(y)

        self.sig_quantile_transformer = self.fit_sig_quantile_transformer(X, y)

        self.set_sample_sexes(X)
        self.set_sample_fusions(X)

        super().fit(X, y)

        self.reset_sample_sexes()
        self.reset_sample_fusions()

        return self

    @property
    def is_fitted(self) -> bool:
        ## Check for an attribute that only exists if the model is fitted
        ## There are many attributes we could check, but one is selected to make the check quick
        return hasattr(self.dna_combined_clf["logistic_regression"], "coef_")

    def _check_is_fitted(self) -> None:
        if not self.is_fitted:
            self.logger.error(self.__class__.__name__ + " is not yet fitted")
            raise Exception


    ## Features ================================
    def get_required_features(self) -> dict[str, pd.Index]:

        self._check_is_fitted()

        required_features = {
            SUB_CLF_NAMES.GEN_POS: self.gen_pos_clf["cos_sim"].profiles_.index,
            SUB_CLF_NAMES.SNV96: self.snv96_clf["logistic_regression"].feature_names_in_,
            SUB_CLF_NAMES.EVENT: self.event_clf["logistic_regression"].feature_names_in_,
            SIG_QUANTILE_TRANSFORMER_NAME: self.sig_quantile_transformer.feature_names_in_,

            SUB_CLF_NAMES.GENE_EXP: self.gene_exp_clf["chi2"].selected_features,
            SUB_CLF_NAMES.ALT_SJ: self.alt_sj_clf["chi2"].selected_features,

            META_CLF_NAMES.DNA_COMBINED: self.dna_combined_clf["logistic_regression"].feature_names_in_,
            META_CLF_NAMES.RNA_COMBINED: self.rna_combined_clf["logistic_regression"].feature_names_in_,
        }

        required_features[META_CLF_NAMES.COMBINED] = (
            required_features[META_CLF_NAMES.DNA_COMBINED].tolist() +
            required_features[META_CLF_NAMES.RNA_COMBINED].tolist()
        )

        required_features = {clf_name: pd.Index(features) for clf_name, features in required_features.items()}

        return required_features

    def _check_features(self, X: pd.DataFrame) -> None:
        handler = MissingFeaturesHandler(X=X, cuppa_classifier=self)
        handler.check_features()

    def fill_missing_cols(
        self,
        X: pd.DataFrame,
        fill_value: int | float = NA_FILL_VALUE,
        verbose: bool = True
    ):

        handler = MissingFeaturesHandler(
            X = X,
            cuppa_classifier=self,
            fill_value = fill_value,
            verbose = verbose
        )

        X_new = handler.fill_missing()

        return X_new


    ## Predicting ================================
    def transform(
        self,
        X: pd.DataFrame,
        y = None,
        until_step: Optional[str] = None,
        keep_steps: Optional[str | list[str]] = None,
        verbose: Optional[bool] = None
    ) -> pd.DataFrame:
        """
        By default, gets the probabilities from the final transformer (i.e. ProbCombiner)

        Can also return probabilities at a specified step (using the `until_step` argument) and/or keep intermediate
        probabilities (using the `keep_steps` argument)

        Parameters
        ----------
        X: pandas DataFrame
            Features (columns) per sample (row)

        y: pandas Series
            Sample labels. Only required when this method is called during fitting (i.e. by fit_transform)

        until_step: str
            Step name to transform until

        keep_steps: str or list[str]
            Which steps to keep? If not None, will return a dict of data frames of transformed features at each step

        verbose: bool
            Show progress messages?
        """

        self._check_is_fitted()
        self._check_features(X)

        self.set_sample_sexes(X)
        self.set_sample_fusions(X)

        X_trans = super().transform(
            X=X,
            y=y,
            until_step=until_step,
            keep_steps=keep_steps,
            verbose=verbose
        )

        self.reset_sample_sexes()
        self.reset_sample_fusions()

        return X_trans

    def predict_proba(
        self,
        X: pd.DataFrame,
        y: pd.Series | None = None,
        verbose: bool = False
    ) -> pd.DataFrame:
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
            Returns the prediction probabilities where:
                Rows are a multi index in the form: (classifier_group, classifier_name, sample_id)
                Columns are the predicted class names

        """

        self._check_is_fitted()
        self._check_features(X)

        keep_steps = [
            LAYER_NAMES.COMBINED,
            LAYER_NAMES.META_CLF,
            LAYER_NAMES.SUB_CLF,
        ]

        probs = self.transform(X, y, keep_steps=keep_steps, verbose=verbose)
        probs = {step: probs[step] for step in keep_steps}

        ## Convert to wide format --------------------------------
        ## Add 'combined__' prefix
        probs[LAYER_NAMES.COMBINED].columns = (
            META_CLF_NAMES.COMBINED +
            DEFAULT_FEATURE_PREFIX_SEPERATOR +
            probs[LAYER_NAMES.COMBINED].columns
        )

        probs = pd.concat(probs.values(), axis=1) ## Shape: n_samples x pred_classes

        ## Make column multi-indexes
        probs.columns = probs.columns.str.split(DEFAULT_FEATURE_PREFIX_SEPERATOR, n=1).map(tuple)
        probs.columns.names = ["clf_name", "pred_class"]

        probs.index.name = "sample_id"

        ## Convert to long format --------------------------------
        probs = probs.stack(level="clf_name") ## Shape: (n_samples, clf_names) x pred_classes

        ## Add classifier group to index
        index = probs.index.to_frame(index=False)

        ## Force row order
        index["sample_id"] = pd.Categorical(index["sample_id"], X.index)
        index["clf_name"] = pd.Categorical(index["clf_name"], CLF_NAMES)
        index["clf_group"] = pd.Categorical(CLF_GROUPS.from_clf_names(index["clf_name"]), CLF_GROUPS.get_all())

        index = index[["sample_id", "clf_group", "clf_name"]]

        probs.index = pd.MultiIndex.from_frame(index)
        probs = probs.sort_index()

        ## Force column order
        probs = probs[self.classes_]

        return probs

    def feat_imp(self) -> pd.DataFrame:
        self._check_is_fitted()
        return FeatureImportance.from_cuppa_classifier(self)

    def feat_contrib(
        self,
        X: pd.DataFrame,
        y: None = None,
        sub_clf_names: str | list[str] | None = None,
        column_type: Literal["features", "classes"] = "features"
    ) -> pd.DataFrame:
        """
        Per sample feature contributions (SHAP values) from the sub-classifiers

        See the below URL for more information on SHAP values:
        https://shap.readthedocs.io/en/latest/example_notebooks/overviews/An%20introduction%20to%20explainable%20AI%20with%20Shapley%20values.html

        For each sample, SHAP values are calculated using the following formula:

            coef[i] * (x[i] - X.mean(axis=0)[i]) for the ith feature

            Where:
            - coef = coefficients of the logistic regression
            - x = feature array of the sample
            - X = feature matrix of the training data where rows are samples and columns are features. `X.mean(axis=0)`
            refers to the mean of each feature across all samples

            Reference: https://shap.readthedocs.io/en/latest/generated/shap.explainers.Linear.html

        Parameters
        ----------
        X: pandas DataFrame
            Values of each feature (columns) for each sample (row)

        y: None
            Not used. Argument only exists for compatibility

        sub_clf_names: str or list of str. Valid values are: 'gen_pos', 'snv96', 'event', 'gene_exp', 'alt_sj'
            If specified, only the feature contributions for the specified sub-classifiers are kept. By default, the
            feature contributions for all sub-classifiers are kept.

        column_type: str, 'features' or 'classes'
            if 'features':
                each column is a multi-index in the form (clf_name, feat_name)
                each row is a multi-index in the form (sample_id, pred_class)

            if 'classes':
                each column is a pred_class
                each row is a multi-index in the form (sample_id, clf_name, feat_name)

        Returns
        -------
        pandas DataFrame
            SHAP values

        """

        self._check_is_fitted()
        self._check_features(X)

        sub_classifiers = self[LAYER_NAMES.SUB_CLF]
        if sub_clf_names is not None:
            super().check_step_names(
                steps = sub_clf_names,
                valid_steps = sub_classifiers.transformer_names,
                logger = self.logger
            )

        feat_contribs = {}
        for name, estimator, columns in sub_classifiers.transformers:
            if sub_clf_names is not None and name not in sub_clf_names:
                continue

            if hasattr(estimator, "feat_contrib"):
                feat_contribs[name] = estimator.feat_contrib(X[columns])

        feat_contribs = pd.concat(feat_contribs, axis=1)
        feat_contribs.columns.names = ["clf_name", "feat_name"]

        if column_type == "features":
            return feat_contribs

        feat_contribs = feat_contribs.unstack(level=0).transpose()
        feat_contribs = feat_contribs.reorder_levels(["sample_id", "clf_name", "feat_name"], axis=0)
        feat_contribs = feat_contribs.loc[X.index]  ## Force original sample order

        return feat_contribs

    def predict(
        self,
        X: pd.DataFrame,
        y: None = None,
        bypass_steps: str | list[str] | None = None,
        verbose: bool = False
    ) -> CuppaPrediction:
        """
        Probabilities and feature contributions/values

        This method outputs a dataframe with the columns:
           sample_id: sample ids
           data_type: 'prob' (probabilities) or 'feat_info' (feature contributions and values)
           clf_group: classifier group
           clf_name: classifier name
           feat_name: feature name. Values are NAs when data_type=='probs'
           feat_value: feature values. Values are NAs when data_type=='probs'
           {pred_class_1}, {pred_class_2}, ...: Each column is a prediction class (i.e. cancer type)

        Parameters
        ----------
        X: pandas DataFrame
           Values of each feature (columns) for each sample (row)

        y: None
           Not used. Argument only exists for compatibility

        verbose: bool
           Show progress messages?

        Returns
        -------
        A pandas DataFrame if export_path is None, else void

        """

        if bypass_steps is None:
            cuppa_classifier = self
        else:
            builder = BypassedClassifierBuilder(self, bypass_steps=bypass_steps, verbose=True)
            cuppa_classifier = builder.build()

        cuppa_classifier._check_is_fitted()

        builder = CuppaPredictionBuilder(cuppa_classifier = cuppa_classifier, X = X, verbose = verbose)
        predictions = builder.build()

        return predictions
