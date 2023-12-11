import pandas as pd
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.misc.mock_data import MockTrainingOutput
from cuppa.classifier.feature_importance import FeatureImportance


class TestFeatureImportance:

    #from cuppa.classifier.test.test_feature_importance import TestFeatureImportance
    #self = TestFeatureImportance

    cuppa_classifier = MockTrainingOutput.cuppa_classifier
    feat_imp = FeatureImportance.from_cuppa_classifier(cuppa_classifier)

    def test_summarize_outputs_expected_columns(self):
        summ = self.feat_imp.summarize()
        expected_cols = pd.Series(["class", "clf_name", "feat_name", "mean", "std", "mean_abs", "rank"])
        assert expected_cols.isin(summ.columns).all()

    def test_get_feat_affixes_outputs_correct_affixes(self):
        feat_names = pd.Series(["gen_pos.Breast", "gen_pos__Breast", "gen_pos.Breast__Other"])
        prefixes, suffixes = self.feat_imp.get_feat_affixes(feat_names, sep = '[.]|__')

        assert prefixes == ["gen_pos", "gen_pos", "gen_pos"]
        assert suffixes == ["Breast", "Breast", "Breast__Other"]