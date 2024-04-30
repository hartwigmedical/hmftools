from tests.mock_data import MockInputData
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, FeatureLoader, FeatureLoader


class TestCuppaFeaturesPaths:

    def test_can_find_files_in_dir_by_pattern(self):
        match_files = CuppaFeaturesPaths.find_files_in_dir_by_pattern(
            directory=MockInputData.dir_new_format,
            pattern="cuppa_data.cohort.snv.*.tsv",
        )
        assert match_files.tolist() == ['cuppa_data.cohort.snv.tsv.gz']


    def test_can_infer_old_features_format_paths_from_dir(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, file_format="old")
        assert True

    def test_can_infer_new_features_format_paths_from_dir(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_new_format, file_format="new")
        assert True


class TestFeatureLoader:

    def test_can_load_single_sample_from_tsv(self):
        loader = FeatureLoader(
            MockInputData.path_tsv_new_format,
            sample_id = "COLO829",
            excl_chroms = ["ChrY", "Y"],
        )

        features = loader.load()
        assert features.shape == (1, 6219)
        assert features.index[0] == "COLO829"

        features_series = features.iloc[0]
        assert features_series["gen_pos.3_26000000"] == 20
        assert features_series["snv96.C>A_ACA"] == 133
        assert features_series["event.tmb.snv_count"] == 37660
        assert features_series["event.sv.SIMPLE_DEL_20KB_1MB"] == 20
        assert features_series.index.str.startswith("event.fusion").sum() == 0
        assert features_series["event.trait.is_male"] == 1
        assert features_series["sig.UV (SBS7)"].round() == 24193

        assert features_series.index.str.startswith("gene_exp").sum() == 0
        assert features_series.index.str.startswith("alt_sj").sum() == 0

    def test_can_load_multi_sample_from_tsvs(self):

        loader = FeatureLoader(MockInputData.dir_new_format)
        features = loader.load()
        assert features.shape == (2, 6225)

        assert features["gen_pos.1_500000"].tolist() == [0, 1]
        assert features["snv96.C>T_TCC"].tolist() == [0, 2]
        assert features["event.tmb.snv_count"].tolist() == [0, 8]
        assert features["event.sv.SIMPLE_DEL_20KB_1MB"].tolist() == [0, 20]
        assert features["event.fusion.TMPRSS2_ERG"].tolist() == [loader.na_fill_value, 1]
        assert features["event.trait.is_male"].tolist() == [0, 1]
        assert features["sig.UV (SBS7)"].tolist() == [0, 6.4]
        assert features["gene_exp.BRAF"].tolist() == [loader.na_fill_value, 3.434]
        assert features["alt_sj.7;140426316;140439612"].tolist() == [loader.na_fill_value, 2]


class TestFeatureLoaderOld:

    def test_can_load_dna_features(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, file_format="old")
        loader = FeatureLoaderOld(paths)
        loader.load_dna_features()
        assert True

    def test_can_load_rna_features(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, file_format="old")
        loader = FeatureLoaderOld(paths)
        loader.load_rna_features()
        assert True

    def test_can_load_from_directory(self):
        paths = CuppaFeaturesPaths.from_dir(MockInputData.dir_old_format, file_format="old")
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
