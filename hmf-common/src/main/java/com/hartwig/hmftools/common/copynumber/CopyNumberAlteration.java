package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.exception.HartwigException;

public enum CopyNumberAlteration {
    GAIN, LOSS, NEUTRAL;

    public static int NORMAL_HUMAN_COPY_NUMBER = 2;
    private static final String IDENTIFIER_ERROR = "Could not parse gain/loss/neutral identifier: %s";

    public static CopyNumberAlteration fromCopyNumber(int value) {
        if (value == NORMAL_HUMAN_COPY_NUMBER) {
            return NEUTRAL;
        }

        return value > NORMAL_HUMAN_COPY_NUMBER ? GAIN : LOSS;
    }

    public static CopyNumberAlteration fromString(String value) throws HartwigException {
        try {
            return CopyNumberAlteration.valueOf(value.toUpperCase());
        } catch (IllegalArgumentException e) {
            throw new HartwigException(String.format(IDENTIFIER_ERROR, value));
        }
    }
}