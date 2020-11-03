package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FusionClassifier {

    private static final Set<String> FUSION_PAIR_KEYWORDS = Sets.newHashSet("Fusion",
            "Disruptive Inframe Deletion",
            "Gene Fusion",
            "fusion",
            "EGFR-KDD",
            "Transcript Regulatory Region Fusion",
            "FGFR3 - BAIAP2L1 Fusion",
            "FLT3-ITD");

    private static final Set<String> PROMISCUOUS_FUSION_KEYWORDS =
            Sets.newHashSet("REARRANGEMENT", "Fusions", "fusion", "rearrange", "Transcript Fusion", "FUSION", "FUSIONS");

    private static final Set<String> INTERNAL_FUSION = Sets.newHashSet("is_deletion", "EGFRvIII", "EGFRvV", "EGFRvII", "ITD");

    private static final Set<String> IGNORED = Sets.newHashSet("3' EXON DELETION", "p61BRAF-V600E");

    private FusionClassifier() {
    }

    public static boolean isFusionPair(@NotNull String featureName, @Nullable String gene, @Nullable String biomarkerType) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene)) {
            return extractFusionEvent(featureName, biomarkerType) == FusionEvent.FUSION_PAIR;
        }

        return false;
    }

    public static boolean isPromiscuousFusion(@NotNull String featureName, @Nullable String gene, @Nullable String biomarkerType) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene)) {
            return extractFusionEvent(featureName, biomarkerType) == FusionEvent.FUSION_PROMISCUOUS;
        }

        return false;
    }

    @Nullable
    public static FusionEvent extractFusionEvent(@NotNull String featureName, @Nullable String biomarkerType) {
        if (IGNORED.contains(featureName)) {
            return null;
        }

        if (INTERNAL_FUSION.contains(featureName)) {
            return FusionEvent.FUSION_PAIR;
        } else if (isTypicalFusionPair(featureName) && !featureName.equals("LOSS-OF-FUNCTION")) {
            return FusionEvent.FUSION_PAIR;
        } else if (FUSION_PAIR_KEYWORDS.contains(featureName)) {
            return FusionEvent.FUSION_PAIR;
        } else if (PROMISCUOUS_FUSION_KEYWORDS.contains(featureName)) {
            if (featureName.contains("-")) {
                return FusionEvent.FUSION_PAIR;
            } else {
                return FusionEvent.FUSION_PROMISCUOUS;
            }
        } else if (biomarkerType != null) {
            if (PROMISCUOUS_FUSION_KEYWORDS.contains(biomarkerType)) {
                if (featureName.contains("-")) {
                    return FusionEvent.FUSION_PAIR;
                } else {
                    return FusionEvent.FUSION_PROMISCUOUS;
                }
            }

            if (FUSION_PAIR_KEYWORDS.contains(biomarkerType)) {
                return FusionEvent.FUSION_PAIR;
            }
        }

        return null;
    }

    private static boolean isTypicalFusionPair(@NotNull String featureName) {
        String trimmedFeature = featureName.trim();
        if (trimmedFeature.contains("-")) {
            String[] parts = trimmedFeature.split("-");
            if (parts.length != 2) {
                return false;
            } else {
                // Assume genes that are fused contain no spaces
                return !parts[0].contains(" ") && !parts[1].contains(" ");
            }
        }
        return false;
    }

    private enum FusionEvent {
        FUSION_PROMISCUOUS,
        FUSION_PAIR
    }
}
