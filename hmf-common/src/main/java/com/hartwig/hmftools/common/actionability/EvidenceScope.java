package com.hartwig.hmftools.common.actionability;

import org.jetbrains.annotations.NotNull;

public enum EvidenceScope {
    BROAD("Broad"),
    SPECIFIC("Specific");

    @NotNull
    private final String readableString;

    EvidenceScope(@NotNull final String readableString) {
        this.readableString = readableString;
    }

    @NotNull
    public String readableString() {
        return readableString;
    }
}
