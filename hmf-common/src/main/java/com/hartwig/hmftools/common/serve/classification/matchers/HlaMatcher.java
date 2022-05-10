package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class HlaMatcher implements EventMatcher{

    @NotNull
    private final Set<String> hlaKeyPhrases;

    HlaMatcher(@NotNull final Set<String> hlaKeyPhrases) {
        this.hlaKeyPhrases = hlaKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : hlaKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }
        return false;
    }
}
