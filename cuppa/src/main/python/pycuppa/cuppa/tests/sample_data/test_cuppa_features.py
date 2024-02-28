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


class TestFeatureLoaderNew:

    #from cuppa.tests.sample_data.test_cuppa_features import TestFeatureLoaderNew
    #self = TestFeatureLoaderNew

    def test_can_load_from_tsv(self):

        loader = FeatureLoaderNew(
            MockInputData.path_tsv_new_format,
            sample_id = "COLO829",
            excl_chroms = ["ChrY", "Y"],
            verbose=True
        )

        features = loader.load()
        assert features.shape[0] == 1
        assert features.shape[1] == 6219
        assert features.index[0] == "COLO829"


        features_series = features.iloc[0]
        assert features_series["gen_pos.3_26000000"] == 20
        assert features_series["snv96.C>A_ACA"] == 133
        assert features_series["event.sv.SIMPLE_DEL_20KB_1MB"] == 20
        assert features_series.index.str.startswith("event.fusion").sum() == 0
        assert features_series["event.trait.is_male"] == 1
        assert features_series["event.tmb.snv_count"] == 37660
        assert features_series["sig.UV (SBS7)"].round() == 24193

        assert features_series.index.str.startswith("gene_exp").sum() == 0
        assert features_series.index.str.startswith("alt_sj").sum() == 0
