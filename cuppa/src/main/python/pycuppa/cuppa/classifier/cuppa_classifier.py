from __future__ import annotations

import os
import shutil
from functools import cached_property
from typing import TYPE_CHECKING, Optional, Literal, Self, Iterable

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import make_column_selector

import cuppa.compose.pipeline
from cuppa.constants import SUB_CLF_NAMES, META_CLF_NAMES, LAYER_NAMES, CLF_GROUPS, DEFAULT_CUPPA_CLASSIFIER_PATH
from cuppa.classifier.classifiers import ClassifierLayers, SubClassifiers, MetaClassifiers
from cuppa.classifier.cuppa_prediction import CuppaPrediction, CuppaPredictionBuilder
from cuppa.classifier.feature_importance import FeatureImportance
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.components.mutational_signatures import SigCohortQuantileTransformer
from cuppa.components.prob_overriders import FusionProbOverrider, SexProbFilter
from cuppa.components.preprocessing import NaRowFilter
from cuppa.compose.column_transformer import DEFAULT_FEATURE_PREFIX_SEPERATOR
from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH, SEX_FEATURE_NAME
from cuppa.logger import LoggerMixin

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
    def to_file(self, path: str, verbose: bool = True) -> None:
        if verbose:
            self.logger.info("Exporting classifier to: " + path)
        joblib.dump(self, filename=path, compress=9)

    @staticmethod
    def _check_is_pickle_file(path: str):
        if not path.endswith(".pickle") and not path.endswith(".pickle.gz"):
            raise ValueError("Path must end with .pickle or .pickle.gz")

    @classmethod
    def from_file(cls, path: str) -> Self:
        cls._check_is_pickle_file(path)
        return joblib.load(path)

    @classmethod
    def copy_model_to_resources(cls, source_path: str, verbose: bool = True) -> None:
        cls._check_is_pickle_file(source_path)

        try:
            cls.from_file(source_path)
        except:
            cls.get_class_logger(cls).info("Could test load classifier. The pickled object is likely not a CuppaClassifier")
            raise IOError

        if verbose:
            cls.get_class_logger(cls).info("Copying classifier to: " + DEFAULT_CUPPA_CLASSIFIER_PATH)

        shutil.copy(source_path, DEFAULT_CUPPA_CLASSIFIER_PATH)

    @classmethod
    def from_resources(cls, verbose: bool = True):
        if verbose:
            cls.get_class_logger(cls).info("Loading classifier from: " + DEFAULT_CUPPA_CLASSIFIER_PATH)
        return cls.from_file(DEFAULT_CUPPA_CLASSIFIER_PATH)

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

    def set_fusion_prob_overrider_bypass(self, bypass_value: bool) -> None:
        if not isinstance(bypass_value, bool):
            self.logger.error("`bypass_value` must be a bool")
            raise TypeError

        self.fusion_prob_overrider.bypass = bypass_value

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

    def enable_prob_cal_bypass(self):
        calibrators = self.prob_calibrators
        for calibrator in calibrators.values():
            calibrator.bypass = True

    def disable_prob_cal_bypass(self):
        calibrators = self.prob_calibrators
        for calibrator in calibrators.values():
            calibrator.bypass = False

    ## Training ================================
    def classes_(self) -> None:
        ## `classes_` is a read-only @property in sklearn.pipeline.Pipeline. Undo this so that it can be overriden
        pass

    def fit_sig_quantile_transformer(self, X: pd.DataFrame, y: pd.Series):
        ## Mutational signature quantiles per class
        ## These are not used for prediction and therefore are generated outside the Pipeline.fit() call
        transformer = SigCohortQuantileTransformer(clip_upper=False)

        transformer = super().fit_cached(
            estimator=transformer,
            X=X[make_column_selector("^sig")],
            y=y,
            cache_path=None if self.fit_cache_dir is None else os.path.join(self.fit_cache_dir, "0-sig_quantile_transformer.pickle.gz"),
            verbose=self.verbose,
            step_name="sig_quantile_transformer",
            logger=self.logger
        )

        return transformer

    def fit(self, X: pd.DataFrame, y: pd.Series) -> Self:

        if self.verbose:
            self.logger.info("Training data: %s samples, %s features, %s classes" % (
                str(X.shape[0]),
                str(X.shape[1]),
                len(y.unique())
            ))

        self.classes_ = np.unique(y)

        self.sig_quantile_transformer = self.fit_sig_quantile_transformer(X, y)

        self.set_sample_sexes(X)
        self.set_sample_fusions(X)

        super().fit(X, y)

        self.reset_sample_sexes()
        self.reset_sample_fusions()

        return self

    @cached_property
    def is_fitted(self) -> bool:
        ## Check for an attribute that only exists if the model is fitted
        ## There are many attributes we could check, but one is selected to make the check quick
        return hasattr(self.dna_combined_clf["logistic_regression"], "coef_")

    def _check_is_fitted(self) -> None:
        if not self.is_fitted:
            self.logger.error(self.__class__.__name__ + " is not yet fitted")
            raise Exception


    ## Features ================================
    @cached_property
    def required_features_by_type(self) -> dict[str, pd.Index]:
        self._check_is_fitted()
        return dict(
            gen_pos = self.gen_pos_clf["cos_sim"].profiles_.index,
            snv96 = pd.Index(self.snv96_clf["logistic_regression"].feature_names_in_),
            event = pd.Index(self.event_clf["logistic_regression"].feature_names_in_),

            sig = self.sig_quantile_transformer.feature_names_in_,

            gene_exp = pd.Index(self.gene_exp_clf["chi2"].selected_features),
            alt_sj = pd.Index(self.alt_sj_clf["chi2"].selected_features)
        )

    @property
    def required_dna_features(self) -> pd.Index:

        features = [
            self.required_features_by_type[feat_type]
            for feat_type in ("gen_pos", "snv96", "event", "sig")
        ]

        return pd.Index(np.concatenate(features))

    @property
    def required_rna_features(self) -> pd.Index:

        features = [
            self.required_features_by_type[feat_type]
            for feat_type in ("gene_exp", "alt_sj")
        ]

        return pd.Index(np.concatenate(features))

    @property
    def required_features(self) -> pd.Index:

        features = [
            self.required_features_by_type[feat_type]
            for feat_type in self.required_features_by_type.keys()
        ]

        return pd.Index(np.concatenate(features))

    @staticmethod
    def _get_feat_type_counts_string(feat_names: pd.Index) -> str:
        feature_types = feat_names.str.extract("(^\w+)", expand=True)

        types, counts = np.unique(feature_types, return_counts=True)
        type_count_string = ", ".join(pd.Series(types) + "=" + pd.Series(counts).astype(str))
        return type_count_string

    def _check_features(self, X: pd.DataFrame) -> None:

        missing = self.required_features[~self.required_features.isin(X.columns)]
        n_missing = len(missing)

        if n_missing == 0:
            return None

        self.logger.error("`X` is missing %i feature columns: %s" % (
            n_missing,
            self._get_feat_type_counts_string(missing)
        ))
        self.logger.error("Please use " + self.__class__.__name__ + ".fill_missing_cols() to ensure `X` has the required columns")
        raise LookupError

    def fill_missing_cols(
        self,
        X: pd.DataFrame,
        fill_value: int | float = 0,
        verbose: bool = True
    ) -> pd.DataFrame:

        ## DNA --------------------------------
        X_dna = X.reindex(columns=self.required_dna_features, fill_value=fill_value)

        ## RNA --------------------------------
        pattern_rna_features = f"^{SUB_CLF_NAMES.GENE_EXP}|{SUB_CLF_NAMES.ALT_SJ}"
        X_rna = X[make_column_selector(pattern_rna_features)]

        if X_rna.shape[1]==0:
            ## If RNA columns are entirely missing, create the RNA matrix of NAs
            X_rna = pd.DataFrame(
                np.full((X.shape[0], len(self.required_rna_features)), np.nan),
                index=X.index,
                columns=self.required_rna_features
            )
        else:
            ## Samples without RNA data either have no RNA columns
            ## We don't want to fill these rows with 0 because this would produce an unwanted probability
            is_missing_rna_data = NaRowFilter.detect_na_rows(X_rna, use_first_col=True)
            X_rna = X_rna\
                .loc[~is_missing_rna_data]\
                .reindex(columns=self.required_rna_features, fill_value=fill_value)

        X_new = pd.concat([X_dna, X_rna], axis=1)
        del X_dna, X_rna

        ## Print info --------------------------------
        removed_features = X.columns[~X.columns.isin(X_new.columns)]
        if verbose and len(removed_features) > 0:
            self.logger.info("Removed %i features that were not required: %s" % (
                len(removed_features),
                self._get_feat_type_counts_string(removed_features)
            ))

        added_features = X_new.columns[~X_new.columns.isin(X.columns)]
        if verbose and len(added_features) > 0:
            self.logger.info("Filled in %i features missing with value=%s: %s" % (
                len(added_features),
                fill_value,
                self._get_feat_type_counts_string(added_features)
            ))

        return X_new


    ## Predicting ================================
    def transform(
        self,
        X: pd.DataFrame,
        y = None,

        until_step: Optional[str] = None,
        keep_steps: Optional[str | list[str]] = None,
        verbose: Optional[bool] = None,

        bypass_prob_cal: bool = False
    ):
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

        bypass_prob_cal: bool
            Skip probability calibration?

        """

        self._check_is_fitted()
        self._check_features(X)

        self.set_sample_sexes(X)
        self.set_sample_fusions(X)
        if bypass_prob_cal:
            self.enable_prob_cal_bypass()

        X_trans = super().transform(
            X=X,
            y=y,
            until_step=until_step,
            keep_steps=keep_steps,
            verbose=verbose
        )

        self.reset_sample_sexes()
        self.reset_sample_fusions()
        if bypass_prob_cal:
            self.disable_prob_cal_bypass()

        return X_trans

    def predict_proba(
        self,
        X: pd.DataFrame,
        y: pd.Series | None = None,
        bypass_prob_cal: bool = False,
        verbose: bool = False
    ) -> pd.DataFrame:
        """
        Parameters
        ----------
        X: pandas DataFrame of shape (n_samples, n_features)
            Where `n_samples` is the number of samples and `n_features` is the number of features

        y: None
            Not used. Argument only exists for compatibility

        bypass_prob_cal: bool
            If True, will return raw probabilities, else calibrated probabilities

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

        probs = self.transform(
            X, y,
            keep_steps=keep_steps,
            bypass_prob_cal=bypass_prob_cal,
            verbose=verbose
        )

        probs = {step: probs[step] for step in keep_steps}

        ## Convert columns to multi-indexes --------------------------------
        ## Add 'combined__' prefix
        probs[LAYER_NAMES.COMBINED].columns = \
            META_CLF_NAMES.COMBINED + \
            DEFAULT_FEATURE_PREFIX_SEPERATOR + \
            probs[LAYER_NAMES.COMBINED].columns.astype(str)

        ## Make columns a list of tuple string pairs (prefix, suffix)
        for probs_i in probs.values():
            probs_i.columns = probs_i.columns.str.split(DEFAULT_FEATURE_PREFIX_SEPERATOR, n=1).map(tuple)

        probs = pd.concat(probs.values(), axis=1)

        probs.columns.names = ["clf_name", "pred_class"]
        probs.index.name = "sample_id"

        ## Long dataframe --------------------------------
        probs_long = probs.stack(level="clf_name") ## Shape: (n_samples, clf_groups, clf_names) x pred_classes

        ## Add classifier group to index
        index = probs_long.index.to_frame(index=False)
        index["clf_group"] = pd.Series(CLF_GROUPS)[index["clf_name"]].values
        index = index[["sample_id", "clf_group", "clf_name"]]

        probs_long.index = pd.MultiIndex.from_frame(index)

        ## Force ordering --------------------------------
        ## Original sample order
        index["sample_id"] = pd.Categorical(index["sample_id"], X.index)

        ## Classifier order
        clf_names = probs.columns.get_level_values("clf_name").unique()
        index["clf_name"] = pd.Categorical(index["clf_name"], clf_names)

        index_reordered = pd.MultiIndex.from_frame(index.sort_values(["sample_id", "clf_name"]))
        probs_long = probs_long.loc[index_reordered]

        ## Class order
        probs_long = probs_long[self.classes_]

        return probs_long

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

        sub_classifiers = self["sub_clfs"]
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
        probs_only: bool = False,
        rm_all_zero_rows: bool = False,
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

        probs_only: bool
            Only return probabilities?

        rm_all_zero_rows: bool
           Per sample, remove features with 0 contribution across all classes? This is intended to reduce the file
           size, especially with many samples

        verbose: bool
           Show progress messages?

        Returns
        -------
        A pandas DataFrame if export_path is None, else void

        """

        self._check_is_fitted()

        builder = CuppaPredictionBuilder(
            cuppa_classifier = self,
            X = X,
            probs_only = probs_only,
            rm_all_zero_rows = rm_all_zero_rows,
            verbose = verbose
        )

        predictions = builder.build()

        return predictions


class MissingFeaturesHandler(LoggerMixin):

    def __init__(
        self,
        X: pd.DataFrame,
        cuppa_classifier: CuppaClassifier | None = None,
        required_features: Iterable | None = None,
        fill_value: int | float = 0,
        verbose: bool = True
    ):
        self.X = X
        self.fill_value = fill_value
        self.verbose = verbose

        self.cuppa_classifier = cuppa_classifier
        if cuppa_classifier is not None:
            self.required_features = self._get_required_features_by_type_from_classifier()
        else:
            self.required_features = pd.Series(required_features)

        self._check_inputs()

    def _check_inputs(self):
        if self.cuppa_classifier is None and self.required_features is None:
            self.logger.error("Either `cuppa_classifier` or `required_features` must be provided")
            raise ValueError

    FEAT_TYPES_DNA = ("gen_pos", "snv96", "event", "sig")
    FEAT_TYPES_RNA = ("gene_exp", "alt_sj")
    FEAT_TYPES = FEAT_TYPES_DNA + FEAT_TYPES_RNA

    @staticmethod
    def _get_required_features_by_type_from_classifier(cuppa_classifier: CuppaClassifier) -> dict[str, pd.Index]:

        cuppa_classifier._check_is_fitted()
        return dict(
            gen_pos=cuppa_classifier.gen_pos_clf["cos_sim"].profiles_.index,
            snv96=pd.Index(cuppa_classifier.snv96_clf["logistic_regression"].feature_names_in_),
            event=pd.Index(cuppa_classifier.event_clf["logistic_regression"].feature_names_in_),

            sig=cuppa_classifier.sig_quantile_transformer.feature_names_in_,

            gene_exp=pd.Index(cuppa_classifier.gene_exp_clf["chi2"].selected_features),
            alt_sj=pd.Index(cuppa_classifier.alt_sj_clf["chi2"].selected_features)
        )

    ## This method is only used for testing
    def _get_required_features_by_type_from_list(self) -> dict[str, pd.Index]:
        required_features = pd.Series(self.required_features)

        d = dict()
        for feat_type in self.FEAT_TYPES:
            d[feat_type] = required_features[required_features.str.startswith(f"{feat_type}.")]

        return d

    @cached_property
    def required_features_by_type(self):

        if self.required_features is not None:
            return self._get_required_features_by_type_from_list()

        if self.cuppa_classifier is not None:
            return self._get_required_features_by_type_from_classifier()

    @property
    def required_dna_features(self) -> pd.Index:

        features = [
            self.required_features_by_type[feat_type]
            for feat_type in ("gen_pos", "snv96", "event", "sig")
        ]

        return pd.Index(np.concatenate(features))

    @property
    def required_rna_features(self) -> pd.Index:

        features = [
            self.required_features_by_type[feat_type]
            for feat_type in ("gene_exp", "alt_sj")
        ]

        return pd.Index(np.concatenate(features))

    @staticmethod
    def _get_feat_type_counts_string(feat_names: pd.Index) -> str:
        feature_types = feat_names.str.extract("(^\w+)", expand=True)

        types, counts = np.unique(feature_types, return_counts=True)
        type_count_string = ", ".join(pd.Series(types) + "=" + pd.Series(counts).astype(str))
        return type_count_string

    def _check_features(self, X: pd.DataFrame) -> None:

        missing = self.required_features[~self.required_features.isin(X.columns)]
        n_missing = len(missing)

        if n_missing == 0:
            return None

        self.logger.error("`X` is missing %i feature columns: %s" % (
            n_missing,
            self._get_feat_type_counts_string(missing)
        ))
        self.logger.error(
            "Please use " + self.__class__.__name__ + ".fill_missing_cols() to ensure `X` has the required columns")
        raise LookupError

    def fill_missing(self) -> pd.DataFrame:

        ## DNA --------------------------------
        X_dna = self.X.reindex(columns=self.required_dna_features, fill_value=self.fill_value)

        ## RNA --------------------------------
        pattern_rna_features = f"^{SUB_CLF_NAMES.GENE_EXP}|{SUB_CLF_NAMES.ALT_SJ}"
        X_rna = self.X[make_column_selector(pattern_rna_features)]

        if X_rna.shape[1] == 0:
            ## If RNA columns are entirely missing, create the RNA matrix of NAs
            X_rna = pd.DataFrame(
                np.full((self.X.shape[0], len(self.required_rna_features)), np.nan),
                index=self.X.index,
                columns=self.required_rna_features
            )
        else:
            ## Samples without RNA data either have no RNA columns
            ## We don't want to fill these rows with 0 because this would produce an unwanted probability
            is_missing_rna_data = NaRowFilter.detect_na_rows(X_rna, use_first_col=True)
            X_rna = X_rna \
                .loc[~is_missing_rna_data] \
                .reindex(columns=self.required_rna_features, fill_value=self.fill_value)

        X_new = pd.concat([X_dna, X_rna], axis=1)
        del X_dna, X_rna

        ## Print info --------------------------------
        removed_features = self.X.columns[~self.X.columns.isin(X_new.columns)]
        if self.verbose and len(removed_features) > 0:
            self.logger.info("Removed %i features that were not required: %s" % (
                len(removed_features),
                self._get_feat_type_counts_string(removed_features)
            ))

        added_features = X_new.columns[~X_new.columns.isin(self.X.columns)]
        if self.verbose and len(added_features) > 0:
            self.logger.info("Filled in %i features missing with value=%s: %s" % (
                len(added_features),
                self.fill_value,
                self._get_feat_type_counts_string(added_features)
            ))

        return X_new