import pandas as pd
import numpy as np

from cuppa.components.preprocessing import Log1pTransformer, MaxScaler, NaRowFilter


class TestLog1pTransformer:
    def test(self):
        X = pd.DataFrame([[2.718-1]])
        X_trans = Log1pTransformer().transform(X)

        assert round(X_trans[0][0], 3) == 1


class TestMaxScaler:
    #from cuppa.components.test.test_preprocessing import TestMaxScaler
    #self=TestMaxScaler

    X_train = pd.DataFrame(
        [[10, 11],
         [-1, -2]],
        columns=["feature_1", "feature_2"],
        index=["sample_1", "sample_2"],
    )

    X_new = pd.DataFrame([
        [100, -100]
    ], columns=["feature_1", "feature_2"])

    def test_with_clipping(self):

        transformer = MaxScaler(clip=True)
        transformer.fit(self.X_train)

        X_new_trans = transformer.transform(self.X_new)

        assert X_new_trans.iloc[0,0] == 1
        assert X_new_trans.iloc[0,1] == 0

    def test_without_clipping(self):

        transformer = MaxScaler(clip=False)
        transformer.fit(self.X_train)

        X_new_trans = transformer.transform(self.X_new)

        assert X_new_trans.iloc[0,0] > 1
        assert X_new_trans.iloc[0,1] < 0


class TestNaRowFilter:

    def test_no_na_rows(self):
        X = pd.DataFrame([
            [0, 0],
            [0, 0]
        ])

        X_trans = NaRowFilter().transform(X)

        assert X.shape == X_trans.shape

    def test_one_na_row(self):
        X = pd.DataFrame([
            [0, 0],
            [0, np.nan]
        ])

        X_trans = NaRowFilter(show_warnings=False).transform(X)

        assert X_trans.shape[0] == 1

    def test_use_first_col_na_first_col(self):
        X = pd.DataFrame([
            [0, 0, 0],
            [np.nan, 0, 0]
        ])

        X_trans = NaRowFilter(show_warnings=False, use_first_col=True).transform(X)

        assert X_trans.shape[0] == 1

    def test_use_first_col_na_not_first_col(self):
        X = pd.DataFrame([
            [0, 0, 0],
            [0, np.nan, 0]
        ])

        X_trans = NaRowFilter(show_warnings=False, use_first_col=True).transform(X)

        assert X_trans.shape[0] == 2

    def test_col_pattern_matching(self):
        X = pd.DataFrame([
            [0, 0, 0],
            [0, np.nan, np.nan]
        ])
        X.columns = ["dna","rna","rna"]

        X_trans = NaRowFilter(show_warnings=False, use_first_col=True, pattern="rna").transform(X)

        assert X_trans.shape[0] == 1
