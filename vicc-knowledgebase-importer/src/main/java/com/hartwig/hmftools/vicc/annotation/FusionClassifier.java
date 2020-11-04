package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FusionClassifier {

    private static final Set<String> FUSION_KEYWORDS = Sets.newHashSet("Disruptive Inframe Deletion",
            "Gene Fusion",
            "Fusion",
            "fusion",
            "FUSION",
            "Fusions",
            "FUSIONS",
            "Transcript Regulatory Region Fusion",
            "Transcript Fusion",
            "REARRANGEMENT",
            "rearrange");

    private static final Set<String> INTERNAL_FUSION_PAIRS =
            Sets.newHashSet("is_deletion", "EGFRvIII", "EGFRvV", "EGFRvII", "EGFR-KDD", "ITD");

    private static final Set<String> FEATURE_NAMES_TO_SKIP = Sets.newHashSet("3' EXON DELETION", "p61BRAF-V600E", "LOSS-OF-FUNCTION");

    private FusionClassifier() {
    }

    public static boolean isFusionPair(@NotNull String featureName, @Nullable String gene, @Nullable String biomarkerType) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene) && !CombinedClassifier.isCombinedEvent(featureName)) {
            return extractFusionEvent(featureName, biomarkerType) == FusionEvent.FUSION_PAIR;
        }

        return false;
    }

    public static boolean isPromiscuousFusion(@NotNull String featureName, @Nullable String gene, @Nullable String biomarkerType) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene) && !CombinedClassifier.isCombinedEvent(featureName)) {
            return extractFusionEvent(featureName, biomarkerType) == FusionEvent.PROMISCUOUS_FUSION;
        }

        return false;
    }

    @Nullable
    public static FusionEvent extractFusionEvent(@NotNull String featureName, @Nullable String biomarkerType) {
        if (FEATURE_NAMES_TO_SKIP.contains(featureName)) {
            return null;
        }

        if (INTERNAL_FUSION_PAIRS.contains(featureName) || isTypicalFusionPair(featureName)) {
            return FusionEvent.FUSION_PAIR;
        } else if (hasFusionKeyword(featureName)) {
            if (featureName.contains("-")) {
                return FusionEvent.FUSION_PAIR;
            } else {
                return FusionEvent.PROMISCUOUS_FUSION;
            }
        }
        //        else if (biomarkerType != null) {
        //            if (FUSION_KEYWORDS.contains(biomarkerType)) {
        //                if (featureName.contains("-")) {
        //                    return FusionEvent.FUSION_PAIR;
        //                } else {
        //                    return FusionEvent.PROMISCUOUS_FUSION;
        //                }
        //            }
        //
        //            if (FUSION_PAIR_KEYWORDS.contains(biomarkerType)) {
        //                return FusionEvent.FUSION_PAIR;
        //            }
        //        }

        return null;
    }

    private static boolean hasFusionKeyword(@NotNull String featureName) {
        for (String keyword : FUSION_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isTypicalFusionPair(@NotNull String featureName) {
        String trimmedFeature = featureName.trim();
        String potentialFusion;
        if (trimmedFeature.contains(" ")) {
            String[] parts = trimmedFeature.split(" ");
            if (!parts[1].equalsIgnoreCase("fusion")) {
                return false;
            }
            potentialFusion = parts[0];
        } else {
            potentialFusion = trimmedFeature;
        }

        if (potentialFusion.contains("-")) {
            String[] parts = potentialFusion.split("-");
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
        PROMISCUOUS_FUSION,
        FUSION_PAIR
    }
}
