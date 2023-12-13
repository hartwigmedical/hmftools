import pandas as pd

from cuppa.misc.mock_data import MockCuppaClassifier, MockTrainingData, MockProbsFromFitTransform
from cuppa.classifier.cuppa_classifier import MissingFeaturesHandler


class TestCuppaClassifier:

    def test_transform_outputs_expected_probs(self):

        classifier = MockCuppaClassifier.cuppa_classifier

        probs = classifier.transform(MockTrainingData.X).round(3)
        probs_expected = MockProbsFromFitTransform.combined.round(3)

        probs = probs.reindex(probs_expected.index).reindex(probs_expected.columns, axis=1)
        assert probs.equals(probs_expected)


class TestMissingFeaturesHandler:

    def test_fill_when_no_rna_features_are_present(self):
        X = pd.DataFrame(
            [[1, 1, 1, 1, 1]],
            columns=["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1", "UNUSED_FEATURE"],
            index=["sample_1"]
        )

        required_dna_features = [
            "gen_pos.feat_1",   "gen_pos.feat_2",
            "snv96.feat_1",     "snv96.feat_2",
            "event.feat_1",     "event.feat_2",
            "sig.feat_1",       "sig.feat_2",
        ]

        required_rna_features = [
            "gene_exp.feat_1",  "gene_exp.feat_2",
            "alt_sj.feat_1",    "alt_sj.feat_2"
        ]

        required_features = pd.Series(required_dna_features + required_rna_features)
        handler = MissingFeaturesHandler(required_features=required_features)

        X_filled = handler.fill_missing_cols(X)

        assert required_features.isin(X_filled.columns).all()
        assert len(handler.removed_features) == 1
        assert len(handler.added_features) == 8

        assert X_filled[required_rna_features].isna().all()
        assert X_filled["gen_pos.feat_2"][0] == 0

    def test_fill_missing_features_from_cuppa_classifier(self):
        classifier = MockCuppaClassifier.cuppa_classifier

        required_features = classifier.required_features ## Mock classifier doesn't use any RNA features

        X = pd.DataFrame(
            [[1, 1, 1, 1]],
            columns=["gen_pos.chr1_0000", "snv96.feat_1", "event.feat_1", "sig.feat_1"],
            index=["sample_1"]
        )

        handler = MissingFeaturesHandler(required_features=required_features)

        X_filled = handler.fill_missing_cols(X)

        assert required_features.isin(X_filled.columns).all()
        assert len(handler.removed_features) == 3

        assert classifier.fill_missing_cols(X).equals(X_filled)
