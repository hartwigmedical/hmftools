package com.hartwig.hmftools.patientreporter.copynumber;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberInterpretation {
    GAIN,
    LOSS,
    NEUTRAL;

    private static final int NORMAL_HUMAN_COPY_NUMBER = 2;

    @NotNull
    public static CopyNumberInterpretation fromCopyNumber(final double value) {
        if (Math.round(value) == NORMAL_HUMAN_COPY_NUMBER) {
            return NEUTRAL;
        }

        return value > NORMAL_HUMAN_COPY_NUMBER ? GAIN : LOSS;
    }
}