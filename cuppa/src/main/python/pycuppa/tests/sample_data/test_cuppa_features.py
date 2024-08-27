import os
import shutil
import tempfile

import pytest

from tests.mock_data import MockInputData
from cuppa.sample_data.cuppa_features import CuppaFeaturesDir, CuppaFeaturesLoader


class TestCuppaFeaturesPaths:

    def test_can_find_files_in_dir_by_pattern(self):
        match_files = CuppaFeaturesDir._find_files_in_dir_by_pattern(
            directory=MockInputData.cohort_dir,
            pattern="cuppa_data.cohort.snv.*.tsv",
        )
        assert match_files.tolist() == ['cuppa_data.cohort.snv.tsv.gz']

    def test_raise_error_when_required_files_missing(self):

        tmp_dir = os.path.join(tempfile.gettempdir(), "CuppaFeaturesPaths_test")

        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

        os.mkdir(tmp_dir)

        with(pytest.raises(FileNotFoundError)):
            CuppaFeaturesDir(tmp_dir).get_file_paths()

        paths = dict(
            ## DNA data is required
            snv = os.path.join(tmp_dir, "cuppa_data.cohort.snv.tsv.gz"),
            sv = os.path.join(tmp_dir, "cuppa_data.cohort.sv.tsv"),
            driver = os.path.join(tmp_dir, "cuppa_data.cohort.driver.tsv"),
            trait = os.path.join(tmp_dir, "cuppa_data.cohort.sample_trait.tsv"),

            ## RNA data is optional
            gene_exp = os.path.join(tmp_dir, "cuppa_data.cohort.gene_exp.tsv.gz"),
            alt_sj = os.path.join(tmp_dir, "cuppa_data.cohort.alt_sj.tsv.gz"),
        )

        ## Check if error is raised if not all DNA data is present
        open(paths["snv"], 'w').close()

        with(pytest.raises(FileNotFoundError)):
            CuppaFeaturesDir(tmp_dir).get_file_paths()

        ## Touch the remaining DNA data and check if CuppaFeaturesPaths succeeds
        open(paths["sv"], 'w').close()
        open(paths["driver"], 'w').close()
        open(paths["trait"], 'w').close()

        CuppaFeaturesDir(tmp_dir).get_file_paths()
        assert True

        shutil.rmtree(tmp_dir)


class TestFeatureLoader:

    def test_can_load_single_sample_from_tsv(self):
        loader = CuppaFeaturesLoader(MockInputData.single_sample_tsv, sample_id ="COLO829")

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

        loader = CuppaFeaturesLoader(MockInputData.cohort_dir)
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
