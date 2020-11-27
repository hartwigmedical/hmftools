package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

class FusionPairAndExonMatcher implements EventMatcher {

    @NotNull
    private final Map<String, Set<String>> fusionPairAndExonsPerGene;

    FusionPairAndExonMatcher(@NotNull final Map<String, Set<String>> fusionPairAndExonsPerGene) {
        this.fusionPairAndExonsPerGene = fusionPairAndExonsPerGene;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        Set<String> entries = fusionPairAndExonsPerGene.get(gene);
        if (entries != null) {
            return entries.contains(event);
        }

        return false;
    }
}
