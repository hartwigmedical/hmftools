package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public enum CopyNumberAlteration {
    GAIN, LOSS, NEUTRAL;

    public static final int NORMAL_HUMAN_COPY_NUMBER = 2;
    private static final String IDENTIFIER_ERROR = "Could not parse gain/loss/neutral identifier: %s";

    public static CopyNumberAlteration fromCopyNumber(final int value) {
        if (value == NORMAL_HUMAN_COPY_NUMBER) {
            return NEUTRAL;
        }

        return value > NORMAL_HUMAN_COPY_NUMBER ? GAIN : LOSS;
    }

    @NotNull
    public static CopyNumberAlteration fromString(@NotNull final String value) throws HartwigException {
        try {
            return CopyNumberAlteration.valueOf(value.toUpperCase());
        } catch (IllegalArgumentException e) {
            throw new HartwigException(String.format(IDENTIFIER_ERROR, value));
        }
    }
}