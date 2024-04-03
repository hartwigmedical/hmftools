import pandas as pd
import pytest

from cuppa.sample_data.sample_metadata import SampleMetadata, TrainingSampleSelector


metadata_df = pd.DataFrame(
    [
        [0, "Breast", "Breast: TNBC", 0],
        [1, "Breast", "Breast: TNBC", 76],
        [2, "Breast", "Breast: TNBC", 151],
        [3, "Breast", "Breast: TNBC", 151],
        [4, "Breast", "Breast: TNBC", 151],

        [5, "Lung", "Lung: Small cell", 151],
        [6, "Lung", "Lung: Small cell", 151],
        [7, "Lung", "Lung: Small cell", 151],
        [8, "Lung", "Lung: Small cell", 151],
        [9, "Lung", "Lung: Small cell", 151],
    ],
    columns=["SampleId", "CancerType", "CancerSubtype", "RnaReadLength"]
)


class TestSampleMetadata:
    # from cuppa.classifier.test.test_data_loader_2 import TestSampleMetadata
    # self=TestSampleMetadata

    def test_init_from_df(self):
        metadata = SampleMetadata.from_data_frame(metadata_df)
        assert isinstance(metadata, SampleMetadata)

    def test_check_required_cols(self):
        df_with_wrong_cols = pd.DataFrame(columns=["col1", "col2"])
        with pytest.raises(KeyError):
            SampleMetadata.from_data_frame(df_with_wrong_cols)


class TestTrainingSampleSelector:

    # from cuppa.sample_data.test.test_sample_metadata import TestTrainingSampleSelector
    # self=TestTrainingSampleSelector

    metadata = SampleMetadata.from_data_frame(metadata_df)

    def test_select_min_rna_samples(self):
        selector = TrainingSampleSelector(
            metadata=self.metadata,
            min_samples_with_rna=5,
            incl_rna_read_lengths=None,
            excl_classes=None
        )

        assert list(selector.selected_samples_rna) == [5,6,7,8,9]

    def test_excl_classes(self):
        selector = TrainingSampleSelector(
            metadata=self.metadata,
            min_samples_with_rna=None,
            incl_rna_read_lengths=None,
            excl_classes="Lung: Small cell"
        )

        assert list(selector.selected_samples_dna) == [0,1,2,3,4]
        assert list(selector.selected_samples_rna) == [0,1,2,3,4]

    def test_incl_rna_read_lengths(self):
        selector = TrainingSampleSelector(
            metadata=self.metadata,
            min_samples_with_rna=None,
            incl_rna_read_lengths=151,
            excl_classes=None
        )

        assert list(selector.selected_samples_rna) == [2,3,4,5,6,7,8,9]

