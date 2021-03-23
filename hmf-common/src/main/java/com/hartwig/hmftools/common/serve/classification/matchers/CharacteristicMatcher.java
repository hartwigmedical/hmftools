package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class CharacteristicMatcher implements EventMatcher {

    @NotNull
    private final Set<String> tumorCharacteristicEvents;

    CharacteristicMatcher(@NotNull final Set<String> tumorCharacteristicEvents) {
        this.tumorCharacteristicEvents = tumorCharacteristicEvents;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        return tumorCharacteristicEvents.contains(event);
    }
}
