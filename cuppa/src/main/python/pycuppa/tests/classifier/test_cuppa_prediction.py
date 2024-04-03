import pandas as pd
import numpy as np

from tests.mock_data import MockCvOutput, MockCuppaClassifier, MockTrainingData
from cuppa.classifier.cuppa_prediction import CuppaPredictionBuilder, CuppaPrediction, CuppaPredSummaryBuilder, CuppaPredSummary


class TestCuppaPredictionBuilder:

    def test_can_build_predictions_from_cuppa_classifier_and_features(self):

        builder = CuppaPredictionBuilder(
            cuppa_classifier=MockCuppaClassifier.cuppa_classifier,
            X=MockTrainingData.X
        )
        predictions = builder.build()

        cancer_types = np.unique(MockTrainingData.y).tolist()
        assert predictions.columns.sort_values().tolist() == cancer_types

        assert list(predictions.index.names) == ['sample_id', 'data_type', 'clf_group', 'clf_name', 'feat_name', 'feat_value']
        assert isinstance(predictions, CuppaPrediction)

    def test_probs_from_builder_equals_probs_from_predict(self):

        cuppa_classifier = MockCuppaClassifier.cuppa_classifier
        X = MockTrainingData.X
        builder = CuppaPredictionBuilder(cuppa_classifier=cuppa_classifier, X=X)

        probs_from_builder = builder.probs.round(3).fillna(0)
        probs_from_predict = cuppa_classifier.predict_proba(X).round(3).fillna(0)

        assert np.all(probs_from_predict.values == probs_from_builder.values)

    def test_corresponding_feature_contrib_and_value_are_on_the_same_row(self):
        cuppa_classifier = MockCuppaClassifier.cuppa_classifier
        X = MockTrainingData.X
        builder = CuppaPredictionBuilder(cuppa_classifier=cuppa_classifier, X=X)

        contribs = builder.feat_contribs

        target_feature = "event.tmb.snv_count"
        target_feature_contribs = contribs[contribs.index.get_level_values("feat_name") == target_feature]
        target_feature_values = target_feature_contribs.loc[X.index].index.get_level_values("feat_value")
        target_feature_values = pd.Series(target_feature_values).astype(int)

        target_feature_values_in_X = X[target_feature]

        assert all(target_feature_values.values == target_feature_values_in_X.values)


class TestCuppaPrediction:

    def test_can_load_from_tsv(self):
        predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
        assert True

    def test_can_cast_from_wide_to_long_dataframe(self):
        predictions_wide = MockCvOutput.predictions
        predictions_long = predictions_wide.get_samples(1).wide_to_long()
        assert True


class TestCuppaPredSummaryBuilder:

    #from cuppa.tests.classifier.test_cuppa_prediction import TestCuppaPredSummaryBuilder
    #self = TestCuppaPredSummaryBuilder

    dummy_predictions = pd.DataFrame.from_records([
        (0, np.nan),
        (0, 0),
        (0, 0),
    ])

    dummy_predictions.columns = ["cancer_type_1", "cancer_type_2"]
    dummy_predictions.index = pd.MultiIndex.from_tuples(
        [
            ("sample_1", "prob", "dna", "dna_combined", np.nan, 0),
            ("sample_1", "feat_contrib", "dna", "event", "event.tmb.snv_count", 0),
            ("sample_1", "sig_quantile", np.nan, np.nan, "SBS7", 0),
        ],
        names=["sample_id", "data_type", "clf_group", "clf_name", "feat_name", "feat_value"]
    )

    dummy_predictions = CuppaPrediction(dummy_predictions)

    def test_can_build_from_dummy_predictions(self):

        builder = CuppaPredSummaryBuilder(predictions=self.dummy_predictions)
        pred_summ = builder.build()

        expected_columns = pd.Series(["sample_id", "clf_name", "pred_class_1", "pred_prob_1"])
        assert expected_columns.isin(pred_summ.columns).all()
        assert isinstance(pred_summ, CuppaPredSummary)

    def test_can_build_from_dummy_predictions_with_actual_classes_provided(self):

        dummy_actual_classes = pd.Series(["cancer_type_1"], index=["sample_1"])

        builder = CuppaPredSummaryBuilder(predictions=self.dummy_predictions, actual_classes=dummy_actual_classes)
        pred_summ = builder.build()

        expected_additional_columns = pd.Series(["actual_class", "is_correct_pred", "which_correct_pred"])
        assert expected_additional_columns.isin(pred_summ.columns).all()


    def test_top_feature_names_and_top_feature_values_dataframes_have_correct_corresponding_values(self):
        predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
        actual_classes = MockTrainingData.y

        ## Build
        builder = CuppaPredSummaryBuilder(predictions=predictions, actual_classes=actual_classes)

        ## Get expected feature values
        target_sample_id = "1_Breast"
        target_pred_class = "Breast"
        target_feat_names = builder.top_features.loc[(target_sample_id, target_pred_class)].values

        expected_feat_values = MockTrainingData.X.loc[target_sample_id, target_feat_names]

        ## Get feature values in pred summary
        actual_feat_values = builder.top_feat_values.loc[(target_sample_id, target_pred_class)]

        assert np.all(expected_feat_values.values == actual_feat_values.values)

    def _test_can_build_from_real_data(self):
        predictions = CuppaPrediction.from_tsv("/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/06-pip_env/cv/report/predictions.tsv.gz")

        metadata = pd.read_csv("/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/cup_ref_sample_data.csv")
        actual_classes = pd.Series(metadata["CancerSubtype"].values, index=metadata["SampleId"])
        actual_classes.index.name = "sample_id"

        builder = CuppaPredSummaryBuilder(predictions, actual_classes=actual_classes, show_extra_info=True, verbose=True)
        self = builder


class TestCuppaPredSummary:

    def test_can_initialize_from_predictions(self):
        predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
        pred_summ = predictions.summarize()
        assert True

    def test_can_load_from_tsv(self):
        pred_summ = CuppaPredSummary.from_tsv(MockCvOutput.path_pred_summ)
        assert True
