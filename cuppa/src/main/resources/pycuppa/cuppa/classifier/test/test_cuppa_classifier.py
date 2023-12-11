from cuppa.constants import DEFAULT_FUSION_OVERRIDES_PATH
from cuppa.misc.mock_data import MockTrainingData, MockProbsFromFitTransform
from cuppa.classifier.cuppa_classifier import CuppaClassifier


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
