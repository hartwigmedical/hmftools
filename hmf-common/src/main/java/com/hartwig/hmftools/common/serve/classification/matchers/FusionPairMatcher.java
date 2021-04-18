package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class FusionPairMatcher implements EventMatcher {

    @NotNull
    private final Set<String> exonicDelDupFusionKeyPhrases;
    @NotNull
    private final Set<String> exonicDelDupFusionEvents;
    @NotNull
    private final Set<String> fusionEventsToSkip;

    public FusionPairMatcher(@NotNull final Set<String> exonicDelDupFusionKeyPhrases, @NotNull final Set<String> exonicDelDupFusionEvents,
            @NotNull final Set<String> fusionEventsToSkip) {
        this.exonicDelDupFusionKeyPhrases = exonicDelDupFusionKeyPhrases;
        this.exonicDelDupFusionEvents = exonicDelDupFusionEvents;
        this.fusionEventsToSkip = fusionEventsToSkip;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (fusionEventsToSkip.contains(event)) {
            return false;
        }

        if (exonicDelDupFusionEvents.contains(event)) {
            return true;
        }

        for (String keyPhrase : exonicDelDupFusionKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

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
            // Assume genes that are fused contain no spaces and do not end with "-"
            return parts.length > 1 && !parts[0].contains(" ") && !parts[1].contains(" ");
        }
        return false;
    }
}
