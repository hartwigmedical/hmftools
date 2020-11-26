package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class PromiscuousFusionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> promiscuousFusionKeywords;
    @NotNull
    private final FusionPairMatcher fusionPairMatcher;

    PromiscuousFusionMatcher(@NotNull final Set<String> promiscuousFusionKeywords,
            @NotNull final FusionPairMatcher fusionPairMatcher) {
        this.promiscuousFusionKeywords = promiscuousFusionKeywords;
        this.fusionPairMatcher = fusionPairMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyword : promiscuousFusionKeywords) {
            if (event.contains(keyword) && !fusionPairMatcher.matches(gene, event)) {
                return true;
            }
        }

        return false;
    }
}
