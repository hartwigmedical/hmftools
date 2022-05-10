package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class CharacteristicMatcher implements EventMatcher {

    @NotNull
    private final Set<String> tumorCharacteristicKeyPhrases;

    CharacteristicMatcher(@NotNull final Set<String> tumorCharacteristicKeyPhrases) {
        this.tumorCharacteristicKeyPhrases = tumorCharacteristicKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : tumorCharacteristicKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }
        return false;
    }
}
