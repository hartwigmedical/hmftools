from cuppa.misc.mock_data import MockTrainingData
from cuppa.components.mutational_signatures import SigCohortQuantileTransformer
from sklearn.compose import make_column_selector

class TestSigCohortQuantileTransformer:
    def test(self):
        X = MockTrainingData.X[make_column_selector("^sig")]
        y = MockTrainingData.y

        transformer = SigCohortQuantileTransformer()
        X_trans = transformer.fit_transform(X, y)

        assert round(X_trans.iloc[0,0], 3) == 0.138