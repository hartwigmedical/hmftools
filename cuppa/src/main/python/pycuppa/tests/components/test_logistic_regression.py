import pandas as pd

from cuppa.components.logistic_regression import LogisticRegression


class TestLogisticRegression:
    X = pd.DataFrame([
        [1, 0, 0],
        [1, 0, 0],
        [1, 0, 0],
        [1, 0, 0],
        [0, 0, 0],

        [0, 1, 0],
        [0, 1, 0],
        [0, 1, 0],
        [0, 1, 0],
        [0, 1, 0],

        [0, 0, 0],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 0],
    ])
    X.columns = ["RUNX1_RUNX1T1", "KIAA1549_BRAF", "TMPRSS2_ERG"]

    y = pd.Series([
        "Acute myeloid leukemia",
        "Acute myeloid leukemia",
        "Acute myeloid leukemia",
        "Acute myeloid leukemia",
        "Acute myeloid leukemia",

        "Liposarcoma",
        "Liposarcoma",
        "Liposarcoma",
        "Liposarcoma",
        "Liposarcoma",

        "Prostate",
        "Prostate",
        "Prostate",
        "Prostate",
        "Prostate",
    ])

    lr = LogisticRegression()
    lr.fit(X, y)

    def test_fit(self):

        coefs = self.lr.feat_imp("coef")

        assert round(coefs.loc["Acute myeloid leukemia","RUNX1_RUNX1T1"], 3) == 1.099
        assert round(coefs.loc["Liposarcoma","KIAA1549_BRAF"], 3) == 1.437
        assert round(coefs.loc["Prostate","TMPRSS2_ERG"], 3) == 0.802

    def test_transform(self):
        X_new = pd.DataFrame([
            [1, 0, 0]
        ])
        X_new.columns = ["RUNX1_RUNX1T1", "KIAA1549_BRAF", "TMPRSS2_ERG"]

        probs = self.lr.predict_proba(X_new)

        assert round(probs.loc[0, "Acute myeloid leukemia"], 3) == 0.725

    def test_feat_contrib(self):
        X_new = pd.DataFrame([
            [1, 0, 0]
        ])
        X_new.columns = ["RUNX1_RUNX1T1", "KIAA1549_BRAF", "TMPRSS2_ERG"]

        contribs = self.lr.feat_contrib(X_new)

        assert round(contribs.loc[(0, "Acute myeloid leukemia"), "RUNX1_RUNX1T1"], 3) == 0.806
