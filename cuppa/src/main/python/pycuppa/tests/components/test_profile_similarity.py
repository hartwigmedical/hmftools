import numpy as np
import pandas as pd

from cuppa.components.profile_similarity import cos_sim, ProfileSimilarityTransformer, NoiseProfileAdder, \
    NonBestSimilarityScaler


class TestCosSim:

    def test_array_vs_array_cos_sim_is_correct(self):
        mut_counts = np.array([1, 0.5, 0.25, 0.1])
        mut_profile = np.array([0.75, 0.65, 0.51, 0.4])
        result = cos_sim(mut_counts, mut_profile)

        assert round(result, 3) == 0.911

    def test_matrix_vs_matrix_cos_sim_values_are_correct(self):

        ## n_samples (2) x n_features (3)
        mut_counts = np.array([
            [105, 8, 21],
            [54, 294, 3]
        ])

        ## n_profiles (4) x n_features (3)
        mut_profiles = np.array([
            [0.9, 0.05, 0.05],
            [0.2, 0.6, 0.2],
            [0.33, 0.33, 0.33],
            [0.1, 0.3, 0.6]
        ])

        ## n_samples (2) x n_profiles (4)
        result = cos_sim(mut_counts, mut_profiles)

        assert result.shape == (2,4)
        assert result.loc[0,0].round(2) == 0.99
        assert result.loc[1,1].round(2) == 0.95


class TestProfileSimilarityTransformer:

    X = pd.DataFrame([
        [10000, 100, 100],
        [10000, 100, 100],
        [10000, 100, 100],

        [10, 20000, 10],
        [10, 20000, 10],
        [10, 20000, 10],
    ])

    y = pd.Series([
        "Breast",
        "Breast",
        "Breast",

        "Lung",
        "Lung",
        "Lung"
    ])

    def test_fit_with_agg_func_sum_gives_expected_profile_matrix(self):
        transformer = ProfileSimilarityTransformer(agg_func="sum")
        transformer.fit(self.X, self.y)

        expected_profiles = pd.DataFrame(
            [[30000, 30],
             [300, 60000],
             [300, 30]
            ]
        )
        expected_profiles.columns = ["Breast", "Lung"]
        assert transformer.profiles_.equals(expected_profiles)

    def test_fit_with_count_ceiling_1000_gives_expected_profile_matrix(self):
        transformer = ProfileSimilarityTransformer(agg_func="sum", count_ceiling=1000)
        transformer.fit(self.X, self.y)

        expected_profiles = pd.DataFrame([
            [2941.2, 1.5],
            [29.4, 2997.0],
            [29.4, 1.5],
        ])
        expected_profiles.columns = ["Breast", "Lung"]

        assert transformer.profiles_.round(1).equals(expected_profiles)


class TestNoiseProfileAdder:
    X: pd.DataFrame = TestProfileSimilarityTransformer.X
    y: pd.Series = TestProfileSimilarityTransformer.y
    transformer = NoiseProfileAdder(agg_func="median", noise_counts=100)

    def test_fit_gives_expected_noise_profile(self):
        self.transformer.fit(self.X, self.y)

        expected_profile = np.array([0.490, 0.504, 0.005])
        assert all(self.transformer.noise_profile.round(3) == expected_profile)

    def test_transform_gives_expected_profile_matrix(self):
        X_trans = self.transformer.fit_transform(self.X, self.y)

        expected_X_trans = pd.DataFrame([
            [10049.0, 150.4, 100.5],
            [10049.0, 150.4, 100.5],
            [10049.0, 150.4, 100.5],
            [59.0, 20050.4, 10.5],
            [59.0, 20050.4, 10.5],
            [59.0, 20050.4, 10.5],
        ])

        assert X_trans.round(1).equals(expected_X_trans)


class TestNonBestSimilarityScaler:
    def test_scaling_gives_correct_cos_sim_values(self):
        cos_sims = pd.DataFrame([
            [0.95, 0.92, 0.87, 0.86, 0.75]
        ])

        transformer = NonBestSimilarityScaler(exponent=5)
        cos_sims_trans = transformer.transform(cos_sims)

        expected_cos_sims_trans = pd.DataFrame([
            [0.95, 0.79, 0.57, 0.54, 0.25]
        ])

        assert cos_sims_trans.round(2).equals(expected_cos_sims_trans)