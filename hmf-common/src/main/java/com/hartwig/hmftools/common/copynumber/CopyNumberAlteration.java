package com.hartwig.hmftools.common.copynumber;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberAlteration {
    GAIN("copy-gain"),
    LOSS("copy-loss"),
    NEUTRAL("none");

    private static final int NORMAL_HUMAN_COPY_NUMBER = 2;

    @NotNull
    private final String description;

    CopyNumberAlteration(@NotNull final String description) {
        this.description = description;
    }

    @NotNull
    public String description() {
        return description;
    }

    @NotNull
    public static CopyNumberAlteration fromCopyNumber(final int value) {
        if (value == NORMAL_HUMAN_COPY_NUMBER) {
            return NEUTRAL;
        }

        return value > NORMAL_HUMAN_COPY_NUMBER ? GAIN : LOSS;
    }
}