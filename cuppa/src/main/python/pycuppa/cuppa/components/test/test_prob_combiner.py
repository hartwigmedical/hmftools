import pandas as pd

from cuppa.misc.mock_data import MockTrainingData, MockTrainingOutput
from cuppa.components.prob_combiner import ProbCombiner


class TestProbCombiner:

    def test_non_zero_probs(self):
        probs = pd.DataFrame(
            [[0.9, 0.8, 0.1, 0.2]],
            columns=[
                "dna_combined__Breast",
                "rna_combined__Breast",
                "dna_combined__Lung",
                "rna_combined__Lung",
            ]
        )

        combiner = ProbCombiner(combine_mode="multiply", prob_floor=0.01)
        probs_combined = combiner.transform(probs).round(3)

        assert probs_combined["Breast"][0] == 0.973

    def test_probs_with_zeros(self):
        probs = pd.DataFrame(
            [[1.0, 0.8, 0.0, 0.2]],
            columns=[
                "dna_combined__Breast",
                "rna_combined__Breast",
                "dna_combined__Lung",
                "rna_combined__Lung",
            ]
        )

        combiner = ProbCombiner(combine_mode="multiply", prob_floor=0.01)
        probs_combined = combiner.transform(probs).round(3)

        assert probs_combined["Breast"][0] == 0.998

    def _probs_from_mock_data(self):
        cuppa_classifier = MockTrainingOutput.cuppa_classifier
        X = MockTrainingData.X
        cuppa_classifier._set_sample_sexes(X)
        dna_rna_combined_probs = cuppa_classifier.transform(X, until_step="meta_clfs")
        return dna_rna_combined_probs


