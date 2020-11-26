package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class PromiscuousFusionMatcher implements EventMatcher {

    private static final Set<String> FUSION_KEYWORDS =
            Sets.newHashSet("Fusion", "fusion", "FUSION", "Fusions", "FUSIONS", "REARRANGEMENT", "rearrange");

    @NotNull
    private final FusionPairMatcher fusionPairMatcher;

    PromiscuousFusionMatcher(@NotNull final FusionPairMatcher fusionPairMatcher) {
        this.fusionPairMatcher = fusionPairMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyword : FUSION_KEYWORDS) {
            if (event.contains(keyword) && !fusionPairMatcher.matches(gene, event)) {
                return true;
            }
        }

        return false;
    }
}
