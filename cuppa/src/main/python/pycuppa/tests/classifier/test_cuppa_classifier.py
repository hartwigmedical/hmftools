import numpy as np
import pandas as pd
import pytest

from cuppa.classifier.cuppa_prediction import CuppaPrediction
from cuppa.constants import PREDICT_NA_FILL_VALUE, CLF_GROUPS
from tests.mock_data import MockTrainingData, MockCuppaClassifier
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_classifier_utils import MissingFeaturesHandler, BypassedClassifierBuilder
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.components.prob_overriders import SexProbFilter, FusionProbOverrider
from cuppa.components.passthrough import PassthroughTransformer


class TestMissingFeaturesHandler:

    def test_missing_rna_features_are_filled_with_small_negative_number(self):
        X = pd.DataFrame(
            [[1, 1, 1, 1]],
            columns=["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1"],
            index=["sample_1"]
        )

        present_features = ["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1"]
        missing_features = ["gen_pos.feat_2", "snv96.feat_2", "event.feat_2", "sig.feat_2"]

        required_features = present_features + missing_features

        fill_value = -1e-8
        handler = MissingFeaturesHandler(X, required_features=required_features, fill_value=fill_value)
        X_filled = handler.fill_missing()

        assert all(X_filled[missing_features].iloc[0] == fill_value)

    def test_missing_rna_features_are_filled_but_only_for_samples_with_rna(self):
        X = pd.DataFrame(
            [
                [1, 1, 1, 1, 1],
                [1, 1, 1, 1, np.nan],
            ],
            columns=["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1", "gene_exp.feat_1",],
            index=["sample_with_rna_data", "sample_without_rna_data"]
        )

        required_dna_features = [
            "gen_pos.feat_1",
            "snv96.feat_1",
            "event.feat_1",
            "sig.feat_1",
        ]

        required_rna_features = [
            "gene_exp.feat_1",  "gene_exp.feat_2",
            "alt_sj.feat_1",    "alt_sj.feat_2"
        ]

        required_features = required_dna_features + required_rna_features

        fill_value = -1e-8
        handler = MissingFeaturesHandler(X, required_features=required_features, fill_value=fill_value)
        X_filled = handler.fill_missing()

        assert X_filled.loc["sample_with_rna_data"].isna().all().__invert__()
        assert X_filled.loc["sample_without_rna_data", required_rna_features].isna().all()

    def test_can_fill_missing_cols_from_classifier(self):
        X = MockTrainingData.X
        y = MockTrainingData.y

        classifier = CuppaClassifier(fusion_overrides_path=None)
        classifier.fit(X, y)

        X_incomplete = X.iloc[:,1:10]

        with pytest.raises(LookupError):
            classifier._check_features(X_incomplete)

        X_filled = classifier.fill_missing_cols(X_incomplete, fill_value=0)
        assert X_filled.shape[1] > X_incomplete.shape[1]


class TestCuppaClassifier:

    classifier = CuppaClassifier(fusion_overrides_path=None)
    classifier.fit(MockTrainingData.X, MockTrainingData.y)

    def test_fit_transform_gives_expected_probabilities(self):
        probs = self.classifier.transform(MockTrainingData.X).round(3)

        probs_expected = MockCuppaClassifier.probs_per_clf["combined"].round(3)
        probs_expected = probs_expected.loc[probs.index, probs.columns]

        assert all(probs == probs_expected)

    @staticmethod
    def _get_feature_group_and_fill_missing(group: str):

        rna_columns = MockTrainingData.X.columns.str.match("^(alt_sj|gene_exp)")

        if group == CLF_GROUPS.DNA:
            features = MockTrainingData.X.loc[:, ~rna_columns]
        elif group == CLF_GROUPS.RNA:
            features = MockTrainingData.X.loc[:, rna_columns]
        else:
            raise ValueError("group must be 'dna' or 'rna'")

        features = TestCuppaClassifier.classifier.fill_missing_cols(features, PREDICT_NA_FILL_VALUE)

        ## Put back sex whch was filtered out by feature selection
        features['event.trait.is_male'] = PREDICT_NA_FILL_VALUE

        return features

    def test_can_predict_with_only_dna_features(self):
        X_dna = self._get_feature_group_and_fill_missing(CLF_GROUPS.DNA)
        prediction = self.classifier.predict(X_dna)
        assert isinstance(prediction, CuppaPrediction)

    def test_can_predict_with_only_rna_features(self):
        X_rna = self._get_feature_group_and_fill_missing(CLF_GROUPS.DNA)
        prediction = self.classifier.predict(X_rna)
        assert isinstance(prediction, CuppaPrediction)


class TestBypassedClassifierBuilder:
    def test_can_replace_transformers_with_passthrough_transformer(self):
        classifier = CuppaClassifier()

        assert isinstance(classifier.sex_filters[0], SexProbFilter)
        assert isinstance(classifier.fusion_prob_overrider, FusionProbOverrider)
        assert isinstance(classifier.prob_calibrators["dna_combined"], RollingAvgCalibration)

        builder = BypassedClassifierBuilder(classifier, bypass_steps="all")

        new_classifier = builder.build()

        assert isinstance(new_classifier.sex_filters[0], PassthroughTransformer)
        assert isinstance(new_classifier.fusion_prob_overrider, PassthroughTransformer)
        assert isinstance(new_classifier.prob_calibrators["dna_combined"], PassthroughTransformer)
