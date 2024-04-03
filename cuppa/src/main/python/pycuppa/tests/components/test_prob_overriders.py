import numpy as np
import pandas as pd

from tests.mock_data import MockTrainingData
from cuppa.components.prob_overriders import SexProbFilter, FusionProbOverrider
from cuppa.classifier.cuppa_classifier import CuppaClassifier


class TestFusionProbOverrider:

    def test_sample_with_two_target_fusions_has_increased_target_probs_after_override(self):

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

        transformer.transform(X=probs, debug=True)

        probs_trans = transformer.transform(X=probs, verbose=True)

        ## Check output
        def _prob_trans_gt_probs(pred_class):
            return probs_trans[pred_class][0] > probs[pred_class][0]

        assert _prob_trans_gt_probs("Head and neck: Adenoid cystic")
        assert _prob_trans_gt_probs("Head and neck: Salivary gland")

    def test_overriden_probs_are_different_to_raw_probs_in_full_cuppa_classifier(self):
        X = MockTrainingData.X
        y = MockTrainingData.y
        fusion_overrides_path = MockTrainingData.path_fusion_overrides

        ## Train classifier
        classifier = CuppaClassifier(fusion_overrides_path=fusion_overrides_path)
        classifier.fit(X, y)

        ## Get probabilities before meta-classifiers
        probs_sub_clf = classifier.transform(X, until_step="sub_clfs")
        probs_sub_clf = probs_sub_clf.loc[:, probs_sub_clf.columns.str.match("gen_pos|snv96|event")]

        ## Get probabilities before 'fusion_overrider' in dna_combined meta-classifier
        dna_combined_clf = classifier.meta_clf_layer["dna_combined"]
        probs_dna_combined = dna_combined_clf.transform(probs_sub_clf, until_step="calibrator")

        overrider = dna_combined_clf["fusion_overrider"]

        ## --------------------------------
        ## Using the mock features as input, the raw and transformed probs are equal because all the non-zero probs
        ## are 1.0. Not useful as a test case...

        overrider.sample_fusions = X.loc[:,X.columns.str.contains("fusion")] ## Set sample_fusions manually. This is normally done by CuppaClassifier.predict()
        transform_output = overrider.transform(probs_dna_combined, debug=True)
        assert transform_output["X"].round(3).equals(transform_output["X_trans"].round(3))

        ## Create dummy test case from probs_dna_combined --------------------------------
        selected_samples = ["97_AML"]

        overrider.sample_fusions = X.loc[selected_samples, X.columns.str.contains("fusion")]

        probs_dna_combined_constructed = pd.DataFrame.from_records(
            [dict(AML=0.5, Breast=0.5,  Lung=0.0,  Melanoma=0.0, Prostate=0.0)],
            index=selected_samples
        )

        transform_output = overrider.transform(probs_dna_combined_constructed, debug=True)

        assert transform_output["X_trans"]["AML"][0] > transform_output["X"]["AML"][0]


class TestSexProbFilter:
    def test_prostate_prob_for_male_sample_is_increased_after_override(self):

        probs = pd.DataFrame(
            [[0.8, 0.2, np.nan]],
            columns=["Prostate", "Ovarian", "Skin"]
        )

        sample_sexes = pd.Series([True]) ## True = male

        probs_trans = SexProbFilter(sample_sexes=sample_sexes).transform(probs)

        assert probs_trans["Prostate"][0] > probs["Prostate"][0]

    def test_ovarian_prob_for_female_sample_is_increased_after_override(self):
        probs = pd.DataFrame(
            [[0.1, 0.9, np.nan]],
            columns=["Prostate", "Ovarian", "Skin"]
        )

        sample_sexes = pd.Series([False]) ## False = female

        probs_trans = SexProbFilter(sample_sexes=sample_sexes).transform(probs)

        assert probs_trans["Ovarian"][0] > probs["Ovarian"][0]