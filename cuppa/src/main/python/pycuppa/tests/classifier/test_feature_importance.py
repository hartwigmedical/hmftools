import numpy as np
import pandas as pd

from tests.mock_data import MockCuppaClassifier, MockTrainingData
from cuppa.classifier.feature_importance import FeatureImportance


class TestFeatureImportance:

    def test_can_initialize_from_cuppa_classifier(self):
        feat_imp = FeatureImportance.from_cuppa_classifier(MockCuppaClassifier.cuppa_classifier)

        cancer_types = np.unique(MockTrainingData.y).tolist()
        assert list(feat_imp.index) == cancer_types
        assert list(feat_imp.columns.names) == ["clf_name", "feat_name"]

    def test_summarize_outputs_expected_columns(self):
        feat_imp = FeatureImportance.from_cuppa_classifier(MockCuppaClassifier.cuppa_classifier)
        summ = feat_imp.summarize()
        expected_cols = pd.Series(["class", "clf_name", "feat_name", "mean", "std", "mean_abs", "rank"])
        assert expected_cols.isin(summ.columns).all()

    def test_get_feat_affixes_outputs_correct_prefixes_and_suffixes(self):

        feat_names = pd.Series(["gen_pos.Breast", "gen_pos__Breast", "gen_pos.Breast__Other"])
        prefixes, suffixes = FeatureImportance.get_feat_affixes(feat_names, sep = '[.]|__')

        assert prefixes == ["gen_pos", "gen_pos", "gen_pos"]
        assert suffixes == ["Breast", "Breast", "Breast__Other"]