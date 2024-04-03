import pandas as pd

from cuppa.components.mutational_signatures import SigCohortQuantileTransformer


class TestSigCohortQuantileTransformer:

    #from cuppa.tests.components.test_mutational_signatures import TestSigCohortQuantileTransformer
    #self = TestSigCohortQuantileTransformer()

    feature_names = ["SBS7_UV", "SBS4_Smoking", "SBS1_Age"]

    sample_ids = [
        "Breast_1", "Breast_2", "Breast_3", "Breast_4", "Breast_5",
        "Lung_1", "Lung_2", "Lung_3", "Lung_4", "Lung_5",
        "Skin_1", "Skin_2", "Skin_3", "Skin_4", "Skin_5",
    ]

    X = pd.DataFrame([

        [2000, 4000, 1000],
        [3000, 2000, 2000],
        [1500, 1000, 1500],
        [3000, 1500, 4000],
        [2000, 3500, 1500],

        [1000, 10000, 4000],
        [2000, 5000, 2000],
        [1500, 15000, 1000],
        [4000, 60000, 1500],
        [1500, 20000, 3500],

        [50000, 1000, 2000],
        [25000, 2000, 3000],
        [40000, 1500, 1500],
        [30000, 3000, 3000],
        [10000, 2000, 2000],
    ], index=sample_ids, columns=feature_names)

    y = pd.Series([
        "Breast",
        "Breast",
        "Breast",
        "Breast",
        "Breast",

        "Lung",
        "Lung",
        "Lung",
        "Lung",
        "Lung",

        "Skin",
        "Skin",
        "Skin",
        "Skin",
        "Skin",
    ], index=sample_ids)

    def test_can_get_clipped_quantiles(self):

        quantile_transformer = SigCohortQuantileTransformer(clip_lower=True, clip_upper=True)
        quantile_transformer.fit(self.X, self.y)

        ## Upper clip
        sample_sig_contribs = pd.DataFrame(
            [dict(SBS7_UV=1_000_000_000, SBS4_Smoking=0, SBS1_Age=0)], ## SBS7 contrib is HIGHER than that found in any cohort defined above
            index=["Skin_new"]
        )

        sample_quantiles = quantile_transformer.transform(sample_sig_contribs)

        expected_sample_quantiles_sbs7 = pd.Series(dict(Breast=1, Lung=1, Skin=1))
        actual_sample_quantiles_sbs7 = sample_quantiles.loc[("Skin_new", "SBS7_UV")].iloc[0]

        assert all(expected_sample_quantiles_sbs7 == actual_sample_quantiles_sbs7)

        ## Lower clip
        sample_sig_contribs = pd.DataFrame(
            [dict(SBS7_UV=1, SBS4_Smoking=0, SBS1_Age=0)],  ## SBS7 contrib is LOWER than that found in any cohort defined above
            index=["Skin_new"]
        )

        sample_quantiles = quantile_transformer.transform(sample_sig_contribs)

        expected_sample_quantiles_sbs7 = pd.Series(dict(Breast=0, Lung=0, Skin=0))
        actual_sample_quantiles_sbs7 = sample_quantiles.loc[("Skin_new", "SBS7_UV")].iloc[0]

        assert all(expected_sample_quantiles_sbs7 == actual_sample_quantiles_sbs7)


    def test_can_unclip_when_sig_contrib_outside_of_cohort_range(self):

        quantile_transformer = SigCohortQuantileTransformer(clip_lower=False, clip_upper=False)
        quantile_transformer.fit(self.X, self.y)

        sample_sig_contribs = pd.DataFrame(
            [dict(SBS7_UV=75000, SBS4_Smoking=1000, SBS1_Age=0)],
            index=["Skin_new"]
        )

        sample_quantiles = quantile_transformer.transform(sample_sig_contribs)

        cohort_quantiles = quantile_transformer.get_quantiles()

        ## Upper unclip
        expected_last_quantile_skin_sbs7_contrib = 50000
        actual_last_quantile_skin_sbs7_contrib = cohort_quantiles.loc["Skin", "SBS7_UV"].iloc[-1]
        assert actual_last_quantile_skin_sbs7_contrib == expected_last_quantile_skin_sbs7_contrib

        expected_quantile_in_sample_sbs7 = 1.5
        actual_quantile_in_sample_skin_sbs7 = sample_quantiles.loc[("Skin_new", "SBS7_UV"), "Skin"].iloc[0]
        assert actual_quantile_in_sample_skin_sbs7 == expected_quantile_in_sample_sbs7

        ## Lower unclip
        expected_first_quantile_lung_sbs4_contrib = 5000
        actual_first_quantile_lung_sbs4_contrib = cohort_quantiles.loc["Lung", "SBS4_Smoking"].iloc[0]
        assert actual_first_quantile_lung_sbs4_contrib == expected_first_quantile_lung_sbs4_contrib

        expected_quantile_in_sample_lung_sbs4 = -5
        actual_quantile_in_sample_lung_sbs4 = sample_quantiles.loc[("Skin_new", "SBS4_Smoking"), "Lung"].iloc[0]
        assert actual_quantile_in_sample_lung_sbs4 == expected_quantile_in_sample_lung_sbs4
