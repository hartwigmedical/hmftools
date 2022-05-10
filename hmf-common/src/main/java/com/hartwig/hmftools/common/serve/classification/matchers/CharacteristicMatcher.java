package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class CharacteristicMatcher implements EventMatcher {

    @NotNull
    private final Set<String> characteristicKeyPhrases;

    CharacteristicMatcher(@NotNull final Set<String> characteristicKeyPhrases) {
        this.characteristicKeyPhrases = characteristicKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : characteristicKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }
        return false;
    }
}
