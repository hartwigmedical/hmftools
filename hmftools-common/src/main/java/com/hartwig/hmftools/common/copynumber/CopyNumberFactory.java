package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public final class CopyNumberFactory {

    private static final String CNV_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int END_FIELD_COLUMN = 2;
    private static final int VALUE_COLUMN = 3;
    private static final int LOSS_GAIN_FIELD_COLUMN = 4;

    private static final String GAIN_IDENTIFIER = "gain";
    private static final String LOSS_IDENTIFIER = "loss";
    private static final String GAIN_LOSS_IDENTIFIER_ERROR = "Could not parse gain/loss identifier: %s";
    private static final String GAIN_LOSS_MISMATCH_WITH_VALUE_ERROR = "gain/loss identifier does not match with value: %s";

    private CopyNumberFactory() {
    }

    @NotNull
    public static CopyNumber fromCNVLine(@NotNull final String cnvLine) throws HartwigException {
        final String[] values = cnvLine.split(CNV_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long start = Long.valueOf(values[START_FIELD_COLUMN].trim()) + 1;
        final long end = Long.valueOf(values[END_FIELD_COLUMN].trim());
        final int value = Integer.valueOf(values[VALUE_COLUMN].trim());
        final String lossOrGain = values[LOSS_GAIN_FIELD_COLUMN].trim();

        if (!(lossOrGain.equals(GAIN_IDENTIFIER) || lossOrGain.equals(LOSS_IDENTIFIER))) {
            throw new HartwigException(String.format(GAIN_LOSS_IDENTIFIER_ERROR, lossOrGain));
        }

        if ((lossOrGain.equals(GAIN_IDENTIFIER) && value <= 2) || lossOrGain.equals(LOSS_IDENTIFIER) && value >= 2) {
            final String identification =
                    chromosome + "[" + Long.toString(start - 1) + " - " + Long.toString(end) + "]";
            throw new HartwigException(String.format(GAIN_LOSS_MISMATCH_WITH_VALUE_ERROR, identification));
        }

        return new CopyNumber(chromosome, start, end, value);
    }
}
