package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

class CombinedMatcher implements EventMatcher {

    @NotNull
    private final Map<String, Set<String>> combinedEventsPerGene;
    @NotNull
    private final HotspotMatcher hotspotMatcher;
    @NotNull
    private final FusionPairMatcher fusionMatcher;
    @NotNull
    private final AmplificationMatcher amplificationMatcher;

    CombinedMatcher(@NotNull final Map<String, Set<String>> combinedEventsPerGene, @NotNull final HotspotMatcher hotspotMatcher,
            @NotNull final FusionPairMatcher fusionMatcher, @NotNull final AmplificationMatcher amplificationMatcher) {
        this.combinedEventsPerGene = combinedEventsPerGene;
        this.hotspotMatcher = hotspotMatcher;
        this.fusionMatcher = fusionMatcher;
        this.amplificationMatcher = amplificationMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        Set<String> entriesForGene = combinedEventsPerGene.get(gene);
        if (entriesForGene != null) {
            if (entriesForGene.contains(event)) {
                return true;
            }
        }

        if ((event.contains(",") || event.contains(";")) && !event.toLowerCase().contains(" or ")) {
            return true;
        } else if (event.contains("+") && !event.toLowerCase().contains("c.") && !event.contains(">")) {
            return true;
        } else if (event.contains("/")) {
            return false;
        } else if (event.trim().contains(" ")) {
            String[] parts = event.trim().replace("  ", " ").split(" ");
            if (fusionMatcher.matches(gene, parts[0])) {
                // Hotspots or amplifications on fusion genes are considered combined.
                return hotspotMatcher.matches(gene, parts[1]) || amplificationMatcher.matches(gene, parts[1]);
            }
        }

        return false;
    }
}
