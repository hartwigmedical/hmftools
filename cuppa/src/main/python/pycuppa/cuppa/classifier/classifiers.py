from __future__ import annotations

import pandas as pd
from sklearn.compose import make_column_selector

from cuppa.constants import SUB_CLF_NAMES, META_CLF_NAMES
from cuppa.components.calibration import RollingAvgCalibration
from cuppa.components.feature_selection import Chi2FeatureSelector, MaxScaledChi2FeatureSelector
from cuppa.components.logistic_regression import LogisticRegression
from cuppa.components.preprocessing import Log1pTransformer, MaxScaler, NaRowFilter
from cuppa.components.prob_combiner import ProbCombiner
from cuppa.components.prob_overriders import FusionProbOverrider, SexProbFilter
from cuppa.components.profile_similarity import ProfileSimilarityTransformer, NonBestSimilarityScaler, NoiseProfileAdder
from cuppa.compose.column_transformer import ColumnTransformer
from cuppa.compose.pipeline import Pipeline

"""
This module contains methods that instantiate components of CUPPA including: 
- the individual classifiers (as sklearn Pipeline-like objects that chain Transformer objects together)
- the sub-classifier and meta-classifier layers (as sklearn ColumnTransformer-like objects)
"""


class SubClassifiers:

    @staticmethod
    def DnaLogisticRegression() -> LogisticRegression:
        return LogisticRegression(
            random_state=0, multi_class="multinomial", class_weight="balanced", solver="saga",
            max_iter=1000, penalty="l1", C=1
        )

    @staticmethod
    def RnaLogisticRegression() -> LogisticRegression:
        return LogisticRegression(
            random_state=0, multi_class="multinomial", class_weight="balanced", solver="saga",
            max_iter=10000, penalty="l1", C=4
        )

    @classmethod
    def GenPosClassifier(cls) -> Pipeline:
        return Pipeline([
            ("add_noise", NoiseProfileAdder(agg_func="sum", count_ceiling=10_000, noise_counts=500)),
            ("cos_sim", ProfileSimilarityTransformer(agg_func="sum", count_ceiling=10_000, feature_prefix=f"{SUB_CLF_NAMES.GEN_POS}.")),
            ("non_best_scale", NonBestSimilarityScaler(exponent=5)),
            ("logistic_regression", cls.DnaLogisticRegression()),
        ])

    @classmethod
    def Snv96Classifier(cls) -> Pipeline:
        return Pipeline([
            ("add_noise", NoiseProfileAdder(agg_func="median", noise_counts=500)),
            ("log1p", Log1pTransformer()),
            ("logistic_regression", cls.DnaLogisticRegression()),
        ])

    ## TODO: add option to whitelist certain features (e.g. event.tmb.snv_count) by Chi2FeatureSelector in EventClassifier
    @classmethod
    def EventClassifier(cls) -> Pipeline:
        return Pipeline([
            ("feature_transform", ColumnTransformer(
                verbose_feature_names_out=False,
                transformers=[
                    ("other",
                     Chi2FeatureSelector(mode="fdr", threshold=0.001),
                     make_column_selector(pattern="^event[.](?:driver|fusion|virus|trait)")),

                    ("sv-tmb",
                     Pipeline([
                         ("log1p", Log1pTransformer()),
                         ("max_scale", MaxScaler(clip=True))
                     ]),
                     make_column_selector(pattern="^event[.](?:sv|tmb)")),
                ]
            )),
            ("logistic_regression", cls.DnaLogisticRegression()),
        ])

    @classmethod
    def GeneExpClassifier(cls) -> Pipeline:
        return Pipeline([
            ("na_filter", NaRowFilter(pattern=f"^{SUB_CLF_NAMES.GENE_EXP}", use_first_col=True, show_warnings=True)),
            ## show_warnings: show removed NA rows once. The remaining `NaRowFilter()` objects will not show warnings
            ("chi2", MaxScaledChi2FeatureSelector(mode="fdr", threshold=0.001, clip=True)),
            ("cos_sim", ProfileSimilarityTransformer(agg_func="mean", feature_prefix=f"{SUB_CLF_NAMES.GENE_EXP}.")),
            ("non_best_scaler", NonBestSimilarityScaler(exponent=5)),
            ("logistic_regression", cls.RnaLogisticRegression()),
        ])

    @classmethod
    def AltSjClassifier(cls) -> Pipeline:
        return Pipeline([
            ("na_filter", NaRowFilter(pattern=f"^{SUB_CLF_NAMES.ALT_SJ}", use_first_col=True, show_warnings=False)),
            # ("logp1", Log1pTransformer()), ## log transform in `DataLoader` beforehand to prevent re-transformation
            ("chi2", MaxScaledChi2FeatureSelector(mode="fdr", threshold=0.001, clip=True)),
            ("cos_sim", ProfileSimilarityTransformer(agg_func="mean", feature_prefix=f"{SUB_CLF_NAMES.ALT_SJ}.")),
            ("non_best_scaler", NonBestSimilarityScaler(exponent=2)),
            ("logistic_regression", cls.RnaLogisticRegression()),
        ])


class MetaClassifiers:

    @staticmethod
    def MetaLogisticRegression() -> LogisticRegression:
        return LogisticRegression(
            random_state=0, multi_class="multinomial", class_weight="balanced", solver="saga",
            max_iter=1000, penalty="l2", C=1
        )

    @staticmethod
    def ProbCalibrator() -> RollingAvgCalibration:
        return RollingAvgCalibration(kernel="gaussian", window_size="variable", min_true_samples=10)

    @classmethod
    def DnaCombinedClassifier(cls, fusion_overrides_path: str | None = None) -> Pipeline:
        fusion_overrides = None
        if fusion_overrides_path is not None:
            fusion_overrides = pd.read_table(fusion_overrides_path)

        return Pipeline([
            ("logistic_regression", cls.MetaLogisticRegression()),
            ("calibrator", cls.ProbCalibrator()),
            ("fusion_overrider", FusionProbOverrider(overrides=fusion_overrides, mask_base_value=0.01)),
            ("sex_filter", SexProbFilter(show_warnings=True)),
        ])

    @classmethod
    def RnaCombinedClassifier(cls) -> Pipeline:

        column_pattern = f"{SUB_CLF_NAMES.GENE_EXP}|{SUB_CLF_NAMES.ALT_SJ}"

        return Pipeline([
            ("na_filter", NaRowFilter(pattern=column_pattern, use_first_col=True, show_warnings=False)),
            ("logistic_regression", cls.MetaLogisticRegression()),
            ("calibrator", cls.ProbCalibrator()),
            ("sex_filter", SexProbFilter(show_warnings=False)),
        ])


class ClassifierLayers:

    @classmethod
    def SubClassifierLayer(cls, n_jobs: int = 1, verbose: bool = False) -> ColumnTransformer:

        return ColumnTransformer(
            n_jobs=n_jobs,
            verbose=verbose,
            verbose_feature_names_out=True,
            transformers=[
                (
                    SUB_CLF_NAMES.GEN_POS,
                    SubClassifiers.GenPosClassifier(),
                    make_column_selector(pattern=f"^{SUB_CLF_NAMES.GEN_POS}")
                ),

                (
                    SUB_CLF_NAMES.SNV96,
                    SubClassifiers.Snv96Classifier(),
                    make_column_selector(pattern=f"^{SUB_CLF_NAMES.SNV96}")
                ),

                (
                    SUB_CLF_NAMES.EVENT,
                    SubClassifiers.EventClassifier(),
                    make_column_selector(pattern=f"^{SUB_CLF_NAMES.EVENT}")
                ),

                (
                    SUB_CLF_NAMES.GENE_EXP,
                    SubClassifiers.GeneExpClassifier(),
                    make_column_selector(pattern=f"^{SUB_CLF_NAMES.GENE_EXP}")
                ),

                (
                    SUB_CLF_NAMES.ALT_SJ,
                    SubClassifiers.AltSjClassifier(),
                    make_column_selector(pattern=f"^{SUB_CLF_NAMES.ALT_SJ}")
                ),
            ]
        )

    @classmethod
    def MetaClassifierLayer(
        cls,
        fusion_overrides_path: str | None = None,
        n_jobs: int = 1,
        verbose: bool = False
    ) -> ColumnTransformer:

        pattern_dna_combined = "^%s|%s|%s" % (
            SUB_CLF_NAMES.GEN_POS,
            SUB_CLF_NAMES.SNV96,
            SUB_CLF_NAMES.EVENT,
        )

        pattern_rna_combined = "^%s|%s" % (
            SUB_CLF_NAMES.GENE_EXP,
            SUB_CLF_NAMES.ALT_SJ,
        )

        return ColumnTransformer(
            n_jobs=n_jobs,
            verbose=verbose,
            verbose_feature_names_out=True,
            transformers=[
                (
                    META_CLF_NAMES.DNA_COMBINED,
                    MetaClassifiers.DnaCombinedClassifier(fusion_overrides_path=fusion_overrides_path),
                    make_column_selector(pattern=pattern_dna_combined)
                ),

                (
                    META_CLF_NAMES.RNA_COMBINED,
                    MetaClassifiers.RnaCombinedClassifier(),
                    make_column_selector(pattern=pattern_rna_combined)
                ),
            ]
        )

    @classmethod
    def CombinedLayer(cls) -> ProbCombiner:
        return ProbCombiner(prob_floor=0.01)
