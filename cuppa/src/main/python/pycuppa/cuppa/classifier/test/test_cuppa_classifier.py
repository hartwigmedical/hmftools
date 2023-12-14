import pandas as pd

from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH
from cuppa.misc.mock_data import MockTrainingData, MockProbsFromFitTransform
from cuppa.classifier.cuppa_classifier import CuppaClassifier, MissingFeaturesHandler


class TestCuppaClassifier:

    def test_fit_transform(self):
        X = MockTrainingData.X
        y = MockTrainingData.y

        classifier = CuppaClassifier(fusion_overrides_path=DEFAULT_FUSION_OVERRIDES_PATH)
        classifier.fit(X, y)

        probs = classifier.transform(X).round(3)
        probs_expected = MockProbsFromFitTransform.combined.round(3)

        probs = probs.reindex(probs_expected.index).reindex(probs_expected.columns, axis=1)
        assert probs.equals(probs_expected)


class TestMissingFeaturesHandler:

    def test_one_sample_without_rna_data(self):
        #X = FeatureLoaderNew(MockInputData.path_tsv_new_format_colo, verbose=True).load()
        X = pd.DataFrame(
            [[1, 1, 1, 1]],
            columns=["gen_pos.feat_1", "snv96.feat_1", "event.feat_1", "sig.feat_1"],
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

        required_features = required_dna_features + required_rna_features

        self = MissingFeaturesHandler(X, required_features=required_features)

        self.fill_missing()
        self.X

        self.required_features


        pass