import pandas as pd
import pytest

from tests.mock_data import MockTrainingData, MockCuppaClassifier
from cuppa.constants import NA_FILL_VALUE
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.classifier.cuppa_classifier_utils import MissingFeaturesHandler, BypassedClassifierBuilder
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.components.prob_overriders import SexProbFilter, FusionProbOverrider
from cuppa.components.passthrough import PassthroughTransformer


class TestCuppaClassifier:

    def test_fit_transform_gives_expected_probabilities(self):
        X = MockTrainingData.X
        y = MockTrainingData.y

        classifier = CuppaClassifier(fusion_overrides_path=None)
        classifier.fit(X, y)

        probs = classifier.transform(X).round(3)

        probs_expected = MockCuppaClassifier.probs_per_clf["combined"].round(3)
        probs_expected = probs_expected.loc[probs.index, probs.columns]

        assert all(probs == probs_expected)


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

        handler = MissingFeaturesHandler(X, required_features=required_features)
        X_filled = handler.fill_missing()

        assert all(X_filled[missing_features].iloc[0] == NA_FILL_VALUE)

    def test_missing_rna_features_are_filled_with_na(self):
        X = pd.DataFrame(
            [[1, 1, 1, 1]],
            columns=["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1"],
            index=["sample_1"]
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

        handler = MissingFeaturesHandler(X, required_features=required_features)
        X_filled = handler.fill_missing()
        assert X_filled[required_rna_features].isna().iloc[0].all()

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
