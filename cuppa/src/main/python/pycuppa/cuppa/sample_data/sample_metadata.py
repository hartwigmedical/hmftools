from functools import cached_property

import pandas as pd
from cuppa.misc.utils import check_required_columns
from cuppa.logger import LoggerMixin


class SampleMetadata(pd.DataFrame, LoggerMixin):

    def __init__(self, df: pd.DataFrame, *args, **kwargs):
        super().__init__(df, *args, **kwargs)

    @property
    def _constructor(self):
        return SampleMetadata

    @staticmethod
    def from_data_frame(df: pd.DataFrame):

        required_columns = ["CancerType", "CancerSubtype", "RnaReadLength"]
        check_required_columns(df, required_columns=required_columns)

        if "SampleId" in df.columns:
            df = df.set_index("SampleId")
            df.index.name = "sample_id"

        return SampleMetadata(df)

    @classmethod
    def from_csv(cls, path: str, *args, **kwargs):
        df = pd.read_csv(path, *args, **kwargs)
        df = cls.from_data_frame(df)
        return SampleMetadata(df)

    @staticmethod
    def from_tsv(path: str):
        return SampleMetadata.from_csv(path, sep='\t')

    def to_tsv(self, path: str, *args, **kwargs):
        self.to_csv(path, sep="\t", *args, **kwargs)


class TrainingSampleSelector(LoggerMixin):

    def __init__(
        self,
        metadata: SampleMetadata,
        min_samples_with_rna: int | None = 5,
        incl_rna_read_lengths: int | list[int] | None = 151,
        excl_classes: str | list[str] | None = ["_Other", "_Unknown"],
        verbose: bool = True
    ):
        self.metadata = metadata.copy()
        self.min_samples_with_rna = min_samples_with_rna
        self.incl_rna_read_lengths = self._get_incl_rna_read_lengths(incl_rna_read_lengths)
        self.excl_classes = pd.Series(excl_classes)
        self.verbose = verbose

    ## Init ================================
    def _get_incl_rna_read_lengths(self, incl_rna_read_lengths):
        if incl_rna_read_lengths is not None:
            return pd.Series(incl_rna_read_lengths)

        return pd.Series(self.metadata["RnaReadLength"].unique())

    def _print_n_excl_samples(self, is_excl_sample: pd.Series, suffix: str = "") -> None:

        if not self.verbose:
            return None

        n_removed_samples = is_excl_sample.sum()
        if n_removed_samples == 0:
            return None

        msg = "Excluded %s samples" % n_removed_samples
        msg += suffix

        self.logger.info(msg)

    ## Inclusions ================================
    def include_classes(self) -> bool | pd.Series:

        if self.excl_classes is None:
            return True

        is_excl_sample = self.metadata["CancerSubtype"].isin(self.excl_classes)

        self._print_n_excl_samples(
            is_excl_sample,
            " with the following types/subtypes: %s " % ", ".join(self.excl_classes)
        )

        is_incl_sample = ~is_excl_sample
        return is_incl_sample

    def include_rna_read_lengths(self) -> bool | pd.Series:

        if self.incl_rna_read_lengths is None:
            return True

        is_incl_sample = self.metadata["RnaReadLength"].isin(self.incl_rna_read_lengths)

        self._print_n_excl_samples(
            ~is_incl_sample,
            " with RNA read lengths other than: %s " % ", ".join(self.incl_rna_read_lengths.astype(str))
        )

        return is_incl_sample

    def include_classes_with_min_rna_samples(self) -> bool | pd.Series:

        if self.min_samples_with_rna is None:
            return True

        subtype_sample_counts = self.metadata.value_counts(["CancerSubtype", "RnaReadLength"]).reset_index()

        excl_classes = subtype_sample_counts.loc[

            (subtype_sample_counts["count"] < self.min_samples_with_rna) &
            subtype_sample_counts["RnaReadLength"].isin(self.incl_rna_read_lengths),

            "CancerSubtype"
        ]
        excl_classes = excl_classes.values
        is_excl_sample = self.metadata["CancerSubtype"].isin(excl_classes)

        self._print_n_excl_samples(
            is_excl_sample,
            " from following types/subtypes with <%s RNA samples: %s " % (
            self.min_samples_with_rna, ", ".join(excl_classes))
        )

        is_incl_sample = ~is_excl_sample
        return is_incl_sample

    ## Main ================================
    @cached_property
    def selected_samples(self) -> dict[str, pd.Series]:

        ## Run inclusions / exclusions
        self.metadata["has_incl_class"] = self.include_classes()
        self.metadata["has_incl_rna_read_length"] = self.include_rna_read_lengths()
        self.metadata["has_class_with_min_rna_samples"] = self.include_classes_with_min_rna_samples()
        # self.metadata

        ## DNA
        self.metadata["is_selected_sample_dna"] = self.metadata["has_incl_class"]

        ## RNA
        self.metadata["is_selected_sample_rna"] = \
            self.metadata["has_incl_class"] & \
            self.metadata["has_incl_rna_read_length"] & \
            self.metadata["has_class_with_min_rna_samples"]

        ##
        selected_samples_dna = self.metadata.index[ self.metadata["is_selected_sample_dna"] ]
        selected_samples_rna = self.metadata.index[ self.metadata["is_selected_sample_rna"] ]

        if self.verbose:
            self.logger.info("Selected %i samples for the DNA classifiers" % len(selected_samples_dna))
            self.logger.info("Selected %i samples for the RNA classifiers" % len(selected_samples_rna))

        return dict(dna=selected_samples_dna, rna=selected_samples_rna)

    @property
    def selected_samples_dna(self) -> pd.Series:
        return self.selected_samples["dna"]

    @property
    def selected_samples_rna(self) -> pd.Series:
        return self.selected_samples["rna"]





