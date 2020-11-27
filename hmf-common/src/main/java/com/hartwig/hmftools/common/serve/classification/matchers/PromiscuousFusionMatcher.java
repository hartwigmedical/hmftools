package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class PromiscuousFusionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> promiscuousFusionKeyPhrases;
    @NotNull
    private final FusionPairMatcher fusionPairMatcher;

    PromiscuousFusionMatcher(@NotNull final Set<String> promiscuousFusionKeyPhrases,
            @NotNull final FusionPairMatcher fusionPairMatcher) {
        this.promiscuousFusionKeyPhrases = promiscuousFusionKeyPhrases;
        this.fusionPairMatcher = fusionPairMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : promiscuousFusionKeyPhrases) {
            if (event.contains(keyPhrase) && !fusionPairMatcher.matches(gene, event)) {
                return true;
            }
        }

        return false;
    }
}
