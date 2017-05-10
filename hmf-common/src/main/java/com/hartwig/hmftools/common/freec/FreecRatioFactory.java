package com.hartwig.hmftools.common.freec;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

public enum FreecRatioFactory {
    ;

    private static final String CNV_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int RATIO_COLUMN = 2;
    private static final int MEDIAN_RATIO_COLUMN = 3;

    @NotNull
    public static List<FreecRatio> loadNormalRatios(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        return loadRatios(FreecFileLoader.normalRatioLines(basePath, sample));
    }

    @NotNull
    public static List<FreecRatio> loadTumorRatios(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        return loadRatios(FreecFileLoader.tumorRatioLines(basePath, sample));
    }

    @NotNull
    static FreecRatio fromRatioLine(@NotNull final String ratioLine) throws HartwigException {
        final String[] values = ratioLine.split(CNV_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long position = Long.valueOf(values[START_FIELD_COLUMN].trim());
        final double ratio = Double.valueOf(values[RATIO_COLUMN].trim());
        final double medianRatio = Double.valueOf(values[MEDIAN_RATIO_COLUMN].trim());

        return ImmutableFreecRatio.of(ratio, medianRatio, chromosome, position);
    }

    @NotNull
    private static List<FreecRatio> loadRatios(@NotNull final List<String> lines) throws IOException, HartwigException {
        final List<FreecRatio> results = Lists.newArrayList();
        for (final String line : lines) {
            FreecRatio ratio = fromRatioLine(line);
            if (ratio.ratio() > -1) {
                results.add(ratio);
            }
        }

        Collections.sort(results);
        return results;
    }
}
