package com.hartwig.hmftools.common.ratio;

import com.hartwig.hmftools.common.exception.HartwigException;
import org.jetbrains.annotations.NotNull;

public enum RatioFactory {
    ;

    private static final String CNV_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int RATIO_COLUMN = 2;
    private static final int MEDIAN_RATIO_COLUMN = 3;

    @NotNull
    public static Ratio fromRatioLine(@NotNull final String ratioLine) throws HartwigException {
        final String[] values = ratioLine.split(CNV_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long position = Long.valueOf(values[START_FIELD_COLUMN].trim()) + 1; // JOBA: Need to verify this!
        final double ratio = Double.valueOf(values[RATIO_COLUMN].trim());
        final double medianRatio = Double.valueOf(values[MEDIAN_RATIO_COLUMN].trim());

        return ImmutableRatio.of(ratio, medianRatio, chromosome, position);
    }
}
