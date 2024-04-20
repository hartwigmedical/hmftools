from __future__ import annotations

import copy
from functools import cached_property
from typing import Iterable, TYPE_CHECKING

import numpy as np
import pandas as pd
from sklearn.compose import make_column_selector

from cuppa.components.passthrough import PassthroughTransformer
from cuppa.components.preprocessing import NaRowFilter
from cuppa.constants import SUB_CLF_NAMES, SIG_QUANTILE_TRANSFORMER_NAME, NA_FILL_VALUE
from cuppa.logger import LoggerMixin

if TYPE_CHECKING:
    from cuppa.classifier.cuppa_classifier import CuppaClassifier


class MissingFeaturesHandler(LoggerMixin):

    def __init__(
        self,
        X: pd.DataFrame,
        cuppa_classifier: CuppaClassifier | None = None,
        required_features: Iterable | None = None,
        fill_value: int | float = NA_FILL_VALUE,
        verbose: bool = True
    ):
        self.X = X
        self.cuppa_classifier = cuppa_classifier
        self.required_features = required_features
        self.fill_value = fill_value
        self.verbose = verbose

        self._check_inputs()

    def _check_inputs(self) -> None:
        if self.cuppa_classifier is None and self.required_features is None:
            self.logger.error("Either `cuppa_classifier` or `required_features` must be provided")
            raise ValueError

    FEAT_TYPES_DNA = SUB_CLF_NAMES.get_dna_clf_names() + [SIG_QUANTILE_TRANSFORMER_NAME]
    FEAT_TYPES_RNA = SUB_CLF_NAMES.get_rna_clf_names()
    FEAT_TYPES = FEAT_TYPES_DNA + FEAT_TYPES_RNA

    def _get_required_features_by_type_from_classifier(self) -> dict[str, pd.Index]:

        required_features = self.cuppa_classifier.get_required_features()
        required_features = {name: required_features[name] for name in self.FEAT_TYPES}

        return required_features

    ## This method is only used for testing
    def _get_required_features_by_type_from_list(self) -> dict[str, pd.Index]:
        required_features = pd.Series(self.required_features)

        d = dict()
        for feat_type in self.FEAT_TYPES:
            d[feat_type] = required_features[required_features.str.startswith(f"{feat_type}.")]

        return d

    @cached_property
    def required_features_by_type(self) -> dict[str, pd.Index]:

        if self.required_features is not None:
            return self._get_required_features_by_type_from_list()

        if self.cuppa_classifier is not None:
            return self._get_required_features_by_type_from_classifier()

    def _get_required_features(self, feat_types: list):
        features = [
            self.required_features_by_type[feat_type]
            for feat_type in feat_types
        ]

        return pd.Index(np.concatenate(features))

    @property
    def required_dna_features(self) -> pd.Index:
        return self._get_required_features(self.FEAT_TYPES_DNA)

    @property
    def required_rna_features(self) -> pd.Index:
        return self._get_required_features(self.FEAT_TYPES_RNA)

    @property
    def required_features_all(self) -> pd.Index:
        return self._get_required_features(self.FEAT_TYPES)

    @staticmethod
    def _get_feat_type_counts_string(feat_names: pd.Index) -> str:
        feature_types = feat_names.str.extract("(^\w+)", expand=True)

        types, counts = np.unique(feature_types, return_counts=True)
        type_count_string = ", ".join(pd.Series(types) + "=" + pd.Series(counts).astype(str))
        return type_count_string

    def check_features(self) -> None:

        missing = self.required_features_all[ ~self.required_features_all.isin(self.X.columns) ]
        n_missing = len(missing)

        if n_missing == 0:
            return None

        self.logger.error("`X` is missing %i feature columns: %s" % (
            n_missing,
            self._get_feat_type_counts_string(missing)
        ))
        self.logger.error("Please use " + self.__class__.__name__ + ".fill_missing_cols() to ensure `X` has the required columns")
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


class BypassedClassifierBuilder(LoggerMixin):
    def __init__(
        self,
        cuppa_classifier: CuppaClassifier,
        bypass_steps: str | list[str],
        inplace: bool = False,
        verbose: bool = True
    ):
        self.cuppa_classifier = cuppa_classifier

        self.bypass_steps = bypass_steps
        self._check_bypass_steps()

        self.inplace = inplace
        self.verbose = verbose

    def _check_bypass_steps(self) -> None:
        valid_bypass_steps = ["calibrator", "fusion_overrider", "sex_filter"]

        if self.bypass_steps == "all":
            self.bypass_steps = valid_bypass_steps

        if ~pd.Series(self.bypass_steps).isin(valid_bypass_steps).all():
            self.logger.error("Valid `bypass_steps` are: " + ", ".join(valid_bypass_steps) + ", or all")
            raise ValueError

    def _bypass_step(self, pipeline, step: str) -> None:
        step_names = [step_i for step_i, transformer in pipeline.steps]
        if step not in step_names:
            self.logger.error("'%s' not found" % step)
            raise LookupError

        step_index = step_names.index(step)
        pipeline.steps[step_index] = (step, PassthroughTransformer())

    def bypass_prob_calibration(self, cuppa_classifier: CuppaClassifier) -> None:
        if self.verbose:
            self.logger.info("Removing probability calibrators")

        self._bypass_step(cuppa_classifier.dna_combined_clf, "calibrator")
        self._bypass_step(cuppa_classifier.rna_combined_clf, "calibrator")

    def bypass_fusion_overrider(self, cuppa_classifier: CuppaClassifier) -> None:
        if self.verbose:
            self.logger.info("Removing fusion overrider")

        self._bypass_step(cuppa_classifier.dna_combined_clf, "fusion_overrider")

    def bypass_sex_filters(self, cuppa_classifier: CuppaClassifier) -> None:
        if self.verbose:
            self.logger.info("Removing sex filter")

        self._bypass_step(cuppa_classifier.dna_combined_clf, "sex_filter")
        self._bypass_step(cuppa_classifier.rna_combined_clf, "sex_filter")

    def build(self):

        if self.inplace:
            if self.verbose:
                self.logger.info("Modifying CuppaClassifier in place")
            cuppa_classifier = self.cuppa_classifier
        else:
            cuppa_classifier = copy.deepcopy(self.cuppa_classifier)

        if "calibrator" in self.bypass_steps:
            self.bypass_prob_calibration(cuppa_classifier)

        if "fusion_overrider" in self.bypass_steps:
            self.bypass_fusion_overrider(cuppa_classifier)

        if "sex_filter" in self.bypass_steps:
            self.bypass_sex_filters(cuppa_classifier)

        return cuppa_classifier
