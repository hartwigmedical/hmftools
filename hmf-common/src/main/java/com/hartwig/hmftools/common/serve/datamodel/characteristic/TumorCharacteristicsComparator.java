package com.hartwig.hmftools.common.serve.datamodel.characteristic;

import org.jetbrains.annotations.NotNull;

public enum TumorCharacteristicsComparator {
    EQUAL_OR_LOWER("<="),
    EQUAL_OR_GREATER(">="),
    LOWER("<"),
    GREATER(">");

    @NotNull
    private final String keyPhrase;

    TumorCharacteristicsComparator(@NotNull final String keyPhrase) {
        this.keyPhrase = keyPhrase;
    }

    @NotNull
    public String keyPhrase() {
        return keyPhrase;
    }
}
