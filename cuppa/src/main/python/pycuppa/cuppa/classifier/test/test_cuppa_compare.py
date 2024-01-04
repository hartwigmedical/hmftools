import pandas as pd

from cuppa.classifier.cuppa_compare import CuppaCompare
from cuppa.classifier.cuppa_prediction import CuppaPredSummary

class TestCuppaCompare:

    def test_align_rows(self):

        df1 = pd.DataFrame([
            [0,"a"],
            [1,"a"],
            [2,"a"],
            [3,"a"]
        ],
            index=[0,1,2,3],
            columns=["id","value"]
        )

        df2 = pd.DataFrame([
            [2,"b"],
            [3,"b"],
            [4,"b"],
        ],
            index=[2,3,4],
            columns=["id", "value"]
        )

        df1_aligned, df2_aligned = CuppaCompare._align_rows(df1, df2, index_columns="id")

        expected_indexes = [0,1,2,3,4]

        assert list(df1_aligned.index) == expected_indexes
        assert list(df2_aligned.index) == expected_indexes