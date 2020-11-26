package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class FusionPairMatcher implements EventMatcher {

    private static final Set<String> EXON_DEL_DUP_FUSION_PAIRS = Sets.newHashSet("EGFRvIII", "EGFRvV", "EGFRvII", "VIII", "EGFR-KDD");

    private static final Set<String> EVENTS_TO_SKIP =
            Sets.newHashSet("AR-V7", "Gain-of-Function", "LOSS-OF-FUNCTION", "LCS6-variant", "DI842-843VM", "FLT3-ITD");

    FusionPairMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (EVENTS_TO_SKIP.contains(event)) {
            return false;
        }

        if (EXON_DEL_DUP_FUSION_PAIRS.contains(event)) {
            return true;
        }  else {
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
    }
}
