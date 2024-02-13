from cuppa.sample_data.cuppa_features import CuppaFeatures
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, FeatureLoaderNew
from cuppa.tests.mock_data import MockInputData


class TestCuppaFeaturesPaths:

    def test_can_infer_old_features_format_paths_from_dir(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, basenames_mode="old")
        assert True

    def test_can_infer_new_features_format_paths_from_dir(self):
        ## TODO: Add test for new cohort level features format
        pass


def assert_features_dataframe_has_all_feature_types(features):
    assert features.columns.str.match("^gen_pos").any()
    assert features.columns.str.match("^snv96").any()
    assert features.columns.str.match("^event[.]sv").any()
    assert features.columns.str.match("^event[.]fusion").any()
    assert features.columns.str.match("^event[.]trait[.]is_male").any()
    assert features.columns.str.match("^event[.]trait[.]whole_genome_duplication").any()
    assert features.columns.str.match("^event[.]tmb").any()
    assert features.columns.str.match("^sig").any()
    assert features.columns.str.match("^gene_exp").any()
    assert features.columns.str.match("^alt_sj").any()


class TestFeatureLoaderOld:

    def test_can_load_dna_features(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, basenames_mode="old")
        loader = FeatureLoaderOld(paths)
        loader.load_dna_features()
        assert True

    def test_can_load_rna_features(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, basenames_mode="old")
        loader = FeatureLoaderOld(paths)
        loader.load_rna_features()
        assert True

    def test_can_load_from_directory(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, basenames_mode="old")
        loader = FeatureLoaderOld(paths)
        features = loader.load_features()
        assert_features_dataframe_has_all_feature_types(features)
        assert True


class TestFeatureLoaderNew:

    #from cuppa.tests.sample_data.test_cuppa_features import TestFeatureLoaderNew
    #self = TestFeatureLoaderNew

    loader = FeatureLoaderNew(MockInputData.path_tsv_new_format_prostate, verbose=True)
    feat_info = loader.parse_feature_names()

    def test_some_signatures_are_excluded(self):
        excluded_features = self.feat_info.query("is_excluded")["key"]
        assert set(excluded_features) == {"SIG_10_POLE", "SIG_11", "SIG_1"}

    def test_expected_renamed_categories_are_present(self):
        categories = self.feat_info["category_renamed"].unique()
        categories = set(categories)

        categories_expected = {
            'gen_pos',
            'snv96',
            # 'driver', ## Test sample is missing drivers
            'event.fusion', 'event.sv', 'event.trait',
            'sig',

            'gene_exp',
            'alt_sj',
        }

        assert categories == categories_expected

    def test_expected_feat_types_are_present(self):
        feat_types = set(self.feat_info["feat_type"].unique())

        feat_types_expected = {
            'gen_pos',
            'snv96',
            # 'driver', ## Test sample is missing drivers
            'event',
            'sig',
            'gene_exp',
            'alt_sj',
        }

        assert feat_types == feat_types_expected

    def test_can_load_from_tsv(self):
        features = FeatureLoaderNew(path=MockInputData.path_tsv_new_format_prostate).load()
        assert_features_dataframe_has_all_feature_types(features)


class TestCuppaFeatures:
    def _test_can_load_from_tsv_files(self):
        ## TODO: Update test when data is available
        directory = "/"
        files = CuppaFeaturesPaths.from_dir(directory, basenames_mode="new")
        features = CuppaFeatures.from_tsv_files(files)
        assert True