import pandas as pd
from cuppa.classifier.cuppa_prediction import CuppaPredictionBuilder, CuppaPrediction, CuppaPredSummaryBuilder, CuppaPredSummary
from cuppa.misc.mock_data import MockCvOutput, MockTrainingOutput, MockTrainingData


class TestCuppaPredictionBuilder:
    #from cuppa.classifier.test.test_cuppa_prediction import TestCuppaPredictionBuilder
    #self=TestCuppaPredictionBuilder

    cuppa_classifier = MockTrainingOutput.cuppa_classifier
    X = MockTrainingData.X
    cuppa_classifier.cv_performance = MockCvOutput.performance

    builder = CuppaPredictionBuilder(cuppa_classifier=cuppa_classifier, X=X)

    def test_probs_builder_equals_probs_predict(self):
        probs_predict = self.cuppa_classifier.predict_proba(self.X)
        probs_builder = self.builder.probs

        assert probs_predict.iloc[0,].round(5).equals( probs_builder.iloc[0,].round(5) )

    def test_feat_contrib_and_values_alignment(self):
        contribs = self.builder.feat_contribs

        sel_feature = "event.tmb.snv_count"

        ## Feat value alignment
        feature_contribs = contribs[contribs.index.get_level_values("feat_name") == sel_feature]

        feature_values = feature_contribs.loc[self.X.index].index.get_level_values("feat_value")
        feature_values = pd.Series(feature_values).astype(int)

        feature_values_X = self.X[sel_feature]

        assert all(feature_values.values == feature_values_X.values)


class TestCuppaPrediction:
    #from cuppa.classifier.test.test_cuppa_prediction import TestCuppaPrediction
    #self = TestCuppaPrediction

    def test_load_from_tsv(self):
        predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
        assert True

    def test_add_cv_performance_yields_the_same_columns(self):
        predictions = MockCvOutput.predictions
        #self=predictions
        perf = MockCvOutput.performance

        predictions_with_perf = predictions.add_cv_performance(perf)

        assert all(predictions.columns == predictions_with_perf.columns)

    def test_wide_to_long_successful(self):
        predictions = MockCvOutput.predictions_for_vis
        predictions.get_samples(1).wide_to_long() #.to_csv("/Users/lnguyen/Desktop/predictions.long.tsv", sep='\t', index=False)
        assert True


class TestCuppaPredSummaryBuilder:
    #from cuppa.classifier.test.test_cuppa_prediction import TestCuppaPredSummaryBuilder
    #self = TestCuppaPredSummaryBuilder

    predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
    actual_classes = MockTrainingData.y
    builder = CuppaPredSummaryBuilder(predictions=predictions, actual_classes=actual_classes)

    def test_alignment_of_top_features_and_top_feat_values(self):
        row = 2
        column = 2

        ##
        top_features = self.builder.top_features
        target_feat_name = top_features.iloc[row, column]
        target_sample_id = top_features.index.get_level_values("sample_id")[row]
        target_pred_class = top_features.index.get_level_values("pred_class")[column]

        ## Get expected feature value
        predictions_index = self.predictions.index.to_frame(index=False)
        expected_feat_value = predictions_index.loc[
            (predictions_index["sample_id"] == target_sample_id) &
            (predictions_index["feat_name"] == target_feat_name),

            "feat_value"
        ]

        ## Get position of target feature
        top_features = self.builder.top_features
        top_features_row = top_features.loc[(target_sample_id, target_pred_class)]

        ## Get target feat value
        top_feat_values = self.builder.top_feat_values
        top_feat_values_row = top_feat_values.loc[(target_sample_id, target_pred_class)]
        target_feat_value = top_feat_values_row[top_features_row == target_feat_name]

        assert expected_feat_value.values[0] == target_feat_value.values[0]


    def test_build_without_feat_info_or_actual_class(self):

        predictions = self.predictions
        predictions = predictions[predictions.index.get_level_values("data_type")=="prob"]

        builder = CuppaPredSummaryBuilder(predictions=predictions, actual_classes=None)
        pred_summ = builder.build()

        assert True

    def _test_build_from_real_data(self):
        predictions = CuppaPrediction.from_tsv("/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/models/Hartwig_PCAWG/29-pre_prod/06-pip_env/cv/report/predictions.tsv.gz")

        metadata = pd.read_csv("/Users/lnguyen/Hartwig/hartwigmedical/analysis/cup/pycuppa/data/features/Hartwig_PCAWG/017/06-NET_as_prefix/tables/cup_ref_sample_data.csv")
        actual_classes = pd.Series(metadata["CancerSubtype"].values, index=metadata["SampleId"])
        actual_classes.index.name = "sample_id"

        builder = CuppaPredSummaryBuilder(predictions, actual_classes=actual_classes, show_extra_info=True, verbose=True)
        self = builder


class TestCuppaPredSummary:

    def test_init_from_predictions(self):
        predictions = CuppaPrediction.from_tsv(MockCvOutput.path_predictions)
        pred_summ = predictions.summarize()
        assert True

    def test_load_from_tsv(self):
        pred_summ = CuppaPredSummary.from_tsv(MockCvOutput.path_pred_summ)
        assert True
