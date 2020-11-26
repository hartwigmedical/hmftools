package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

class FusionPairAndExonRangeMatcher implements EventMatcher {

    @NotNull
    private final Map<String, Set<String>> fusionPairAndExonRangesPerGene;

    FusionPairAndExonRangeMatcher(@NotNull final Map<String, Set<String>> fusionPairAndExonRangesPerGene) {
        this.fusionPairAndExonRangesPerGene = fusionPairAndExonRangesPerGene;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        Set<String> entries = fusionPairAndExonRangesPerGene.get(gene);
        if (entries != null) {
            return entries.contains(event);
        }

        return false;
    }
}
