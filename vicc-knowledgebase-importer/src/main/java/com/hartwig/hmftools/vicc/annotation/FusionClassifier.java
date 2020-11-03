package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FusionClassifier {

    //TODO: check if more EXON_DEL_DUP fusions need to be added
    private static final Set<String> FUSION_PAIR_AND_EXON_RANGE =
            Sets.newHashSet("KIT EXON 11 MUTATION", "KIT Exon 11 mutations", "KIT Exon 11 deletions", "MET EXON 14 SKIPPING MUTATION");

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

    private static final Set<String> IGNORE = Sets.newHashSet("3' EXON DELETION");

    private FusionClassifier() {
    }

    public static boolean isFusionPairAndGeneRangeExon(@Nullable String featureDescription) {
        return FUSION_PAIR_AND_EXON_RANGE.contains(featureDescription);
    }

    public static boolean isFusionPair(@NotNull String featureName, @Nullable String biomarkerType) {
        return extractFusionEvent(featureName, biomarkerType) == FusionEvent.FUSION_PAIR;
    }

    public static boolean isPromiscuousFusion(@NotNull String featureName, @Nullable String biomarkerType) {
        return extractFusionEvent(featureName, biomarkerType) == FusionEvent.FUSION_PROMISCUOUS;
    }

    @Nullable
    private static FusionEvent extractFusionEvent(@NotNull String featureName, @Nullable String biomarkerType) {
        if (!IGNORE.contains(featureName)) { // Extract internal fusion
            if (INTERNAL_FUSION.contains(featureName) || featureName.contains("[a-zA-Z]+")) {
                return FusionEvent.FUSION_PAIR;
            } else if (featureName.contains("-") && !featureName.equals("LOSS-OF-FUNCTION")) {
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
        }

        return null;
    }

    private enum FusionEvent {
        FUSION_PROMISCUOUS,
        FUSION_PAIR
    }
}
