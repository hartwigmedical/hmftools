import pandas as pd

from cuppa.components.preprocessing import MaxScaler, Log1pTransformer
from cuppa.compose.column_transformer import ColumnTransformer


class TestColumnTransformer:

    def test_column_transformer_outputs_correct_values(self):
        X = pd.DataFrame([
            [0, 0, 0, 2, 2],
            [1, 1, 1, 2, 2],
            [9, 9, 9, 2, 2],
        ])
        X.columns = ["feat_1","feat_2","feat_3","feat_4","feat_5"]
        X.index = ["sample_1", "sample_2", "sample_3"]

        column_transformer = ColumnTransformer(transformers=[
            ("group_1", MaxScaler(), ["feat_1","feat_2","feat_3"]),
            ("group_2", Log1pTransformer(), ["feat_4","feat_5"]),
        ])

        column_transformer.fit(X)
        X_trans = column_transformer.transform(X, verbose_feature_names_out=False)

        assert X_trans.loc["sample_3","feat_1"] == 1.0
        assert X_trans.loc["sample_2","feat_1"].round(3) == 0.111

        assert X_trans.loc["sample_1","feat_4"].round(3) == 1.099

