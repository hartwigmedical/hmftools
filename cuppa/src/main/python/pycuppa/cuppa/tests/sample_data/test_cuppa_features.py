from cuppa.sample_data.cuppa_features import CuppaFeatures
from cuppa.classifier.cuppa_classifier import CuppaClassifier
from cuppa.misc.mock_data import MockInputData
from cuppa.sample_data.cuppa_features import CuppaFeaturesPaths, FeatureLoaderOld, FeatureLoaderNew


class TestCuppaFeaturesPaths:

    def _test_from_dir(self):
        #self = CuppaDataFiles(dict())

        # directory = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables"
        # files = CuppaFeaturesPaths.from_dir(directory)

        directory = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/output/cuppa_features"
        files = CuppaFeaturesPaths.from_dir(directory, basenames_mode="old")
        assert True

class TestCuppaFeatures:

    def _test_from_csv_files_old_format(self):
        directory = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables"
        paths = CuppaFeaturesPaths.from_dir(directory, basenames_mode="old")

        loader = FeatureLoaderOld(paths)
        X = loader.load_features()

    def _test_from_tsv_new_format_predict(self):
        model_path = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/04-NET_as_prefix/cuppa_classifier.pickle.gz"
        cuppa_classifier = CuppaClassifier.from_file(model_path)

        features_path = MockInputData.path_tsv_new_format_colo

        loader = FeatureLoaderNew(features_path)
        features = loader.load()
        features = cuppa_classifier.fill_missing_cols(features)

        predictions = cuppa_classifier.predict(features, verbose=True)
        pred_summ = predictions.summarize()

        assert pred_summ.loc[pred_summ["clf_group"]=="combined", "pred_class_1"].values == "Prostate"

    def _test_from_tsv_files(self):
        #cls = CuppaFeatures(pd.DataFrame())
        directory = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/output/cuppa_features"
        files = CuppaFeaturesPaths.from_dir(directory, basenames_mode="new")
        X = CuppaFeatures.from_tsv_files(files)

        assert True


class TestFeatureLoaderOld:
    input_dir = "/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/"
    paths_features = CuppaFeaturesPaths.from_dir(input_dir, basenames_mode="old")
    loader = FeatureLoaderOld(paths_features)

    def _test_load_dna_features(self):
        X = self.loader.load_dna_features()
        assert True

    def _test_load_rna_features(self):
        X = self.loader.load_rna_features()
        assert True


class TestFeatureLoaderNew:

    #loader = FeatureLoaderNew(MockInputData.path_tsv_new_format_colo, verbose=True)
    loader = FeatureLoaderNew(MockInputData.path_tsv_new_format_prostate, verbose=True)
    feat_info = loader.parse_feature_names()

    def test_exclude_features(self):
        assert self.feat_info.query("is_excluded").shape[0] == 3

    def test_categories_correct(self):
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

    def test_feat_types_correct(self):
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

    def test_load_successful(self):
        X = self.loader.load()
        assert True
