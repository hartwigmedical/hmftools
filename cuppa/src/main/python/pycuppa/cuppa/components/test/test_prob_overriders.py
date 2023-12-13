import numpy as np
import pandas as pd
from cuppa.components.prob_overriders import SexProbFilter, FusionProbOverrider


class TestFusionProbOverrider:

    def test_sample_with_two_target_fusions(self):
        # from cuppa.components.postprocessing import FusionProbOverrider
        # import pandas as pd

        ## Inputs
        probs = pd.DataFrame(
            [[0.35, 0.30, 0.25, 0.10]],
            columns=["Head and neck: Adenoid cystic", "Head and neck: Salivary gland", "Breast", "Prostate"]
        )

        sample_fusions = pd.DataFrame(
            [[1.0, 1.0]],
            columns=["event.fusion.MYB_NFIB", "event.fusion.MYBL1_NFIB"]
        )

        overrides = pd.DataFrame(
            [
                ["event.fusion.", "MYB_NFIB", "Head and neck: Adenoid cystic"],
                ["event.fusion.", "MYBL1_NFIB", "Head and neck: Salivary gland"]
            ],
            columns=["feat_prefix","feat_basename","target_class"]
        )

        ## Transform
        transformer = FusionProbOverrider(
            sample_fusions=sample_fusions,
            overrides=overrides
        )

        probs_trans = transformer.transform(probs)

        ## Check output
        def _prob_trans_gt_probs(pred_class):
            return probs_trans[pred_class][0] > probs[pred_class][0]

        assert _prob_trans_gt_probs("Head and neck: Adenoid cystic")
        assert _prob_trans_gt_probs("Head and neck: Salivary gland")


class TestSexProbFilter:
    def test_male_sample(self):

        probs = pd.DataFrame(
            [[0.8, 0.2, np.nan]],
            columns=["Prostate", "Ovarian", "Skin"]
        )

        sample_sexes = pd.Series([True]) ## True = male

        probs_trans = SexProbFilter(sample_sexes=sample_sexes).transform(probs)

        assert probs_trans["Prostate"][0] > probs["Prostate"][0]

    def test_female_sample(self):
        probs = pd.DataFrame(
            [[0.1, 0.9, np.nan]],
            columns=["Prostate", "Ovarian", "Skin"]
        )

        sample_sexes = pd.Series([False]) ## False = female

        probs_trans = SexProbFilter(sample_sexes=sample_sexes).transform(probs)

        assert probs_trans["Ovarian"][0] > probs["Ovarian"][0]