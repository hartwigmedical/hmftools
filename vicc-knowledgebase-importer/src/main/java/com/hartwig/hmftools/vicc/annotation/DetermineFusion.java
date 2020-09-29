package com.hartwig.hmftools.vicc.annotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class DetermineFusion {

    private DetermineFusion() {
    }

    public static boolean isFusion(@NotNull String feature, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @NotNull String proteinAnnotation) {
        FusionEvent eventKeyFusion = extractKeyFusion(feature, biomarkerType, provenanceRule, proteinAnnotation);
        if (eventKeyFusion.equals(FusionEvent.FUSION_PAIR)) {
            return true;
        } else {
            return false;
        }
    }

    public static boolean isFusionPromiscuous(@NotNull String feature, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @NotNull String proteinAnnotation) {
        FusionEvent eventKeyFusion = extractKeyFusion(feature, biomarkerType, provenanceRule, proteinAnnotation);
        if (eventKeyFusion.equals(FusionEvent.FUSION_PROMISCUOUS)) {
            return true;
        } else {
            return false;
        }
    }

    @NotNull
    private static FusionEvent extractKeyFusion(@NotNull String feature, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @NotNull String proteinAnnotation) {
        if (!FeatureTypeExtractor.IGNORE.contains(feature)) { // Extract internal fusion
            if (FeatureTypeExtractor.INTERNAL_FUSION.contains(proteinAnnotation)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("DELETION") && FeatureTypeExtractor.INTERNAL_FUSION.contains(biomarkerType)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("DELETION") && FeatureTypeExtractor.INTERNAL_FUSION.contains(provenanceRule)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("deletion")) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("deletion/insertion")) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("insertions/deletions")) {
                return FusionEvent.FUSION_PAIR;
            } else if (FeatureTypeExtractor.SEARCH_FUSION_PAIRS.contains(proteinAnnotation)) {
                return FusionEvent.FUSION_PAIR;
            } else if (FeatureTypeExtractor.SEARCH_FUSION_PAIRS.contains(feature)) {
                return FusionEvent.FUSION_PAIR;
            } else if (FeatureTypeExtractor.SEARCH_FUSION_PROMISCUOUS.contains(proteinAnnotation)) {
                if (feature.contains("-")) {
                    return FusionEvent.FUSION_PAIR;
                } else {
                    return FusionEvent.FUSION_PROMISCUOUS;
                }
            } else if (biomarkerType != null) {
                if (FeatureTypeExtractor.SEARCH_FUSION_PROMISCUOUS.contains(biomarkerType)) {
                    if (feature.contains("-")) {
                        return FusionEvent.FUSION_PAIR;
                    } else {
                        return FusionEvent.FUSION_PROMISCUOUS;
                    }
                }
                if (FeatureTypeExtractor.SEARCH_FUSION_PAIRS.contains(biomarkerType)) {
                    return FusionEvent.FUSION_PAIR;
                }
            } else if (feature.toLowerCase().contains("exon") && feature.toLowerCase().contains("deletion") && feature.contains("-")
                    || feature.contains("&")) {// Extract internal fusion
                return FusionEvent.FUSION_PAIR;
            } else {
                return FusionEvent.UNKNOWN;
            }
        }

        //TODO: check why this is needed??
        return FusionEvent.UNKNOWN;
    }
}
