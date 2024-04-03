import pandas as pd

from cuppa.components.prob_combiner import ProbCombiner


class TestProbCombiner:

    def test_can_multiply_probs_without_zeroes(self):
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

    def test_can_multiply_probs_with_zeroes(self):
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

