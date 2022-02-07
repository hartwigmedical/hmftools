package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class WildTypeMatcher implements EventMatcher{
    @NotNull
    private final Set<String> geneWildTypesKeyPhrases;

    public WildTypeMatcher(@NotNull final Set<String> geneWildTypesKeyPhrases) {
        this.geneWildTypesKeyPhrases = geneWildTypesKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : geneWildTypesKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
