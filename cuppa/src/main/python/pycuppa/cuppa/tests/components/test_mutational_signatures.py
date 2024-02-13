import numpy as np
import pandas as pd
from cuppa.components.mutational_signatures import SigCohortQuantileTransformer


class TestSigCohortQuantileTransformer:
    def test_quantiles_are_correct_for_selected_sample(self):

        X = pd.DataFrame([
            [2000, 4000, 1000],
            [3000, 2000, 2000],
            [1500, 1000, 1500],
            [3000, 1500, 4000],
            [2000, 3500, 1500],

            [1000, 10000, 4000],
            [2000, 5000,  2000],
            [1500, 15000, 1000],
            [4000, 60000, 1500],
            [1500, 20000, 3500],

            [50000, 1000, 2000],
            [25000, 2000, 3000],
            [40000, 1500, 1500],
            [30000, 3000, 3000],
            [10000, 2000, 2000],
        ])

        feature_names = ["SBS7_UV", "SBS4_Smoking", "SBS1_Age"]

        sample_ids = [
            "Breast_1", "Breast_2", "Breast_3", "Breast_4", "Breast_5",
            "Lung_1", "Lung_2", "Lung_3", "Lung_4", "Lung_5",
            "Skin_1", "Skin_2", "Skin_3", "Skin_4", "Skin_5",
        ]

        X.columns = feature_names
        X.index = sample_ids

        y = pd.Series([
            "Breast", "Breast", "Breast", "Breast", "Breast",
            "Lung", "Lung", "Lung", "Lung", "Lung",
            "Skin", "Skin", "Skin", "Skin", "Skin",

        ], index=sample_ids)

        transformer = SigCohortQuantileTransformer()
        transformer.fit(X, y)

        sigs_one_sample = pd.DataFrame(
            [[35000, 0, 0]],
            columns=feature_names,
            index=["Skin_6"]
        )

        quantiles_one_sample = transformer.transform(sigs_one_sample)\
            .round(3)\
            .loc["Skin_6", "SBS7_UV"]\
            .values[0].tolist()

        assert np.all(quantiles_one_sample == [1.0, 1.0, 0.625])
