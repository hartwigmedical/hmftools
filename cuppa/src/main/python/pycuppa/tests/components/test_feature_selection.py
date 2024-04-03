import pandas as pd
from sklearn.compose import make_column_selector

from tests.mock_data import MockTrainingData
from cuppa.components.feature_selection import Chi2FeatureSelector


class TestChi2FeatureSelector:

    X = pd.DataFrame([
        [10, 0, 0, 0],
        [10, 1, 0, 0],
        [10, 0, 0, 0],

        [0, 0, 5, 0],
        [0, 1, 5, 0],
        [0, 0, 5, 0],

        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
    ], columns=["feat_1", "feat_2", "feat_3", "feat_4"])

    y = [
        "class_1",
        "class_1",
        "class_1",

        "class_2",
        "class_2",
        "class_2",

        "class_3",
        "class_3",
        "class_3",
    ]

    selector = Chi2FeatureSelector(mode="qvalue", threshold=0.1)
    selector.fit(X, y)

    def test_feat_1_results_are_correct(self):

        results = self.selector.test_results

        chi2_stat = results.loc["feat_1", "stat"]
        assert chi2_stat == 60

        pvalue = results.loc["feat_1", "pvalue"]
        assert "{0:.3g}".format(pvalue) == '9.36e-14'

    def test_rank_is_correct(self):
        assert list(self.selector.test_results["rank"]) == [1, 3, 2, 4]

    def test_transform_returns_correct_features(self):
        X_trans = self.selector.transform(self.X)
        assert list(X_trans.columns) == ["feat_1", "feat_3"]

    def _test_with_event_features(self):
        X = MockTrainingData.X[make_column_selector(pattern="^event[.](?:driver|fusion|virus|trait)")]
        y = MockTrainingData.y
        selector = Chi2FeatureSelector(mode="fdr", threshold=0.001)
        selector.fit(X, y)



