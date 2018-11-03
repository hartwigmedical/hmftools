package com.hartwig.hmftools.common.actionability.cnv;

import org.jetbrains.annotations.NotNull;

enum CopyNumberType {
    AMPLIFICATION("Amplification"),
    DELETION("Deletion");

    @NotNull
    static CopyNumberType fromString(@NotNull String type) {
        if (type.equalsIgnoreCase(AMPLIFICATION.readableString())) {
            return AMPLIFICATION;
        } else if (type.equalsIgnoreCase(DELETION.readableString())) {
            return DELETION;
        }

        throw new IllegalArgumentException("Could not resolve copy number type: " + type);
    }

    @NotNull
    private final String readableString;

    CopyNumberType(@NotNull final String readableString) {
        this.readableString = readableString;
    }

    @NotNull
    public String readableString() {
        return readableString;
    }
}
