package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FusionClassifier {

    private static final Set<String> FUSION_KEYWORDS =
            Sets.newHashSet("Fusion", "fusion", "FUSION", "Fusions", "FUSIONS", "REARRANGEMENT", "rearrange");

    private static final Set<String> EXON_DEL_DUP_FUSION_PAIRS = Sets.newHashSet("EGFRvIII", "EGFRvV", "EGFRvII", "VIII", "EGFR-KDD");

    private static final Set<String> EVENTS_TO_SKIP =
            Sets.newHashSet("AR-V7", "Gain-of-Function", "LOSS-OF-FUNCTION", "LCS6-variant", "DI842-843VM", "FLT3-ITD");

    private FusionClassifier() {
    }

    public static boolean isFusionPair(@NotNull String gene, @NotNull String event) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(gene, event) && !CombinedClassifier.isCombinedEvent(gene, event)) {
            return extractFusionEvent(event) == FusionEvent.FUSION_PAIR;
        }

        return false;
    }

    public static boolean isPromiscuousFusion(@NotNull String gene, @NotNull String event) {
        if (!CombinedClassifier.isFusionPairAndGeneRangeExon(gene, event) && !CombinedClassifier.isCombinedEvent(gene, event)) {
            return extractFusionEvent(event) == FusionEvent.PROMISCUOUS_FUSION;
        }

        return false;
    }

    @Nullable
    public static FusionEvent extractFusionEvent(@NotNull String event) {
        if (EVENTS_TO_SKIP.contains(event)) {
            return null;
        }

        if (EXON_DEL_DUP_FUSION_PAIRS.contains(event) || isTypicalFusionPair(event)) {
            return FusionEvent.FUSION_PAIR;
        } else if (hasFusionKeyword(event)) {
            return FusionEvent.PROMISCUOUS_FUSION;
        }

        return null;
    }

    private static boolean hasFusionKeyword(@NotNull String event) {
        for (String keyword : FUSION_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isTypicalFusionPair(@NotNull String event) {
        String trimmedEvent = event.trim();
        String potentialFusion;
        if (trimmedEvent.contains(" ")) {
            String[] parts = trimmedEvent.split(" ");
            if (!parts[1].equalsIgnoreCase("fusion")) {
                return false;
            }
            potentialFusion = parts[0];
        } else {
            potentialFusion = trimmedEvent;
        }

        if (potentialFusion.contains("-")) {
            String[] parts = potentialFusion.split("-");
            // Assume genes that are fused contain no spaces
            return !parts[0].contains(" ") && !parts[1].contains(" ");
        }
        return false;
    }

    private enum FusionEvent {
        PROMISCUOUS_FUSION,
        FUSION_PAIR
    }
}
