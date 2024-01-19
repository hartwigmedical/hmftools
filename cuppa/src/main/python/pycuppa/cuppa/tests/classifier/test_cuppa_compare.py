import numpy as np
import pandas as pd

from cuppa.classifier.cuppa_compare import CuppaCompare
from cuppa.classifier.cuppa_prediction import CuppaPredSummary


class TestCuppaCompare:

    #from cuppa.tests.classifier.test_cuppa_compare import TestCuppaCompare
    #self=TestCuppaCompare

    pred_summ_old = pd.DataFrame.from_records([
        dict(
            sample_id="sample_1",
            actual_class="cancer_type_1",
            clf_group="dna",
            clf_name="dna_combined",
            is_correct_pred=True,
            which_correct_pred=1,
            pred_class_1="cancer_type_1",
            pred_class_2="cancer_type_2",
            pred_prob_1=0.99,
            pred_prob_2=0.01,
        )
    ])

    pred_summ_new = pd.DataFrame.from_records([
        dict(
            sample_id="sample_1",
            actual_class="cancer_type_1",
            clf_group="dna",
            clf_name="dna_combined",
            is_correct_pred=False,
            which_correct_pred=2,
            pred_class_1="cancer_type_2",
            pred_class_2="cancer_type_1",
            pred_prob_1=0.65,
            pred_prob_2=0.35,
        ),

        dict(
            sample_id="sample_2",
            actual_class="cancer_type_3",
            clf_group="dna",
            clf_name="dna_combined",
            is_correct_pred=False,
            which_correct_pred=1,
            pred_class_1="cancer_type_3",
            pred_class_2="cancer_type_4",
            pred_prob_1=0.90,
            pred_prob_2=0.10,
        )
    ])

    pred_summ_old = CuppaPredSummary.from_data_frame(pred_summ_old)
    pred_summ_new = CuppaPredSummary.from_data_frame(pred_summ_new)

    comparer = CuppaCompare(
        pred_summ_old=pred_summ_old,
        pred_summ_new=pred_summ_new
    )

    def test_compare_predictions_has_expected_columns(self):
        expected_columns = pd.MultiIndex.from_tuples([
            (    'info',       'sample_id'),
            (    'info',        'clf_name'),
            (    'info',    'correct_type'),

            ('is_equal',    'actual_class'),
            ('is_equal',    'pred_class_1'),
            ('is_equal',     'pred_prob_1'),

            (     'new',    'actual_class'),
            (     'new',    'pred_class_1'),
            (     'new',     'pred_prob_1'),
            (     'new', 'is_correct_pred'),

            (     'old',    'actual_class'),
            (     'old',    'pred_class_1'),
            (     'old',     'pred_prob_1'),
            (     'old', 'is_correct_pred'),
        ])

        actual_columns = self.comparer.compare_predictions().columns

        assert expected_columns.equals(actual_columns)

    def test_can_compare_predictions_from_dummy_data(self):

        comparison = self.comparer.compare_predictions()

        correct_types = comparison[("info", "correct_type")]
        assert correct_types[0] == "old_only"
        assert np.isnan(correct_types[1])

        def common_column_values_are_equal(df1: pd.DataFrame, df2: pd.DataFrame) -> bool:
            common_columns = df1.columns.intersection(df2.columns)

            return np.all(
                df1[common_columns].values ==
                df2[common_columns].values
            )

        assert common_column_values_are_equal(
            comparison[comparison[("info", "sample_id")] == "sample_1"]["new"],
            self.pred_summ_new.query("sample_id=='sample_1'")
        )

        assert common_column_values_are_equal(
            comparison[comparison[("info", "sample_id")] == "sample_2"]["new"],
            self.pred_summ_new.query("sample_id=='sample_2'")
        )

        assert common_column_values_are_equal(
            comparison[comparison[("info", "sample_id")] == "sample_1"]["old"],
            self.pred_summ_old.query("sample_id=='sample_1'")
        )

        assert np.all(
            comparison[comparison[("info", "sample_id")] == "sample_2"]["old"].isna().values
        )

    def test_can_compare_performance_has_expected_columns(self):
        expected_columns = pd.MultiIndex.from_tuples([
            ('info', 'class'),
            ('info', 'clf_name'),

            ('diff', 'n_total'),
            ('diff', 'n_predicted'),
            ('diff', 'n_correct'),
            ('diff', 'recall'),
            ('diff', 'precision'),

            ('new', 'n_total'),
            ('new', 'n_predicted'),
            ('new', 'n_correct'),
            ('new', 'recall'),
            ('new', 'precision'),

            ('old', 'n_total'),
            ('old', 'n_predicted'),
            ('old', 'n_correct'),
            ('old', 'recall'),
            ('old', 'precision'),
        ])

        actual_columns = self.comparer.compare_performance().columns

        assert expected_columns.equals(actual_columns)

    def test_align_rows_method_gives_two_dataframes_with_the_same_indexes(self):

        df1 = pd.DataFrame([
            [0,"a"],
            [1,"a"],
            [2,"a"],
            [3,"a"]
        ],
            index=[0,1,2,3],
            columns=["id","value"]
        )

        df2 = pd.DataFrame([
            [2,"b"],
            [3,"b"],
            [4,"b"],
        ],
            index=[2,3,4],
            columns=["id", "value"]
        )

        df1_aligned, df2_aligned = CuppaCompare._align_rows(df1, df2, index_columns="id")

        expected_indexes = [0,1,2,3,4]

        assert list(df1_aligned.index) == expected_indexes
        assert list(df2_aligned.index) == expected_indexes