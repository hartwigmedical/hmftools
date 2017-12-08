package com.hartwig.hmftools.common.copynumber;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberAlteration {
    GAIN,
    LOSS,
    NEUTRAL;

    private static final int NORMAL_HUMAN_COPY_NUMBER = 2;

    @NotNull
    public static CopyNumberAlteration fromCopyNumber(final int value) {
        if (value == NORMAL_HUMAN_COPY_NUMBER) {
            return NEUTRAL;
        }

        return value > NORMAL_HUMAN_COPY_NUMBER ? GAIN : LOSS;
    }
}