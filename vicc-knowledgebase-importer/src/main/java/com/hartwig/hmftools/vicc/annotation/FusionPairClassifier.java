package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

class FusionPairClassifier implements EventClassifier {

    private static final Set<String> EXON_DEL_DUP_FUSION_PAIRS = Sets.newHashSet("EGFRvIII", "EGFRvV", "EGFRvII", "VIII", "EGFR-KDD");

    private static final Set<String> EVENTS_TO_SKIP =
            Sets.newHashSet("AR-V7", "Gain-of-Function", "LOSS-OF-FUNCTION", "LCS6-variant", "DI842-843VM", "FLT3-ITD");

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new FusionPairClassifier());
    }

    private FusionPairClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (EVENTS_TO_SKIP.contains(event)) {
            return false;
        }

        return EXON_DEL_DUP_FUSION_PAIRS.contains(event) || isTypicalFusionPair(event);
    }

    public static boolean isTypicalFusionPair(@NotNull String event) {
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
