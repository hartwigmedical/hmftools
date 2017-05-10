package com.hartwig.hmftools.common.freec;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.common.exception.HartwigException;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.util.List;

public enum FreecCopyNumberFactory {
    ;

    private static final String CNV_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int END_FIELD_COLUMN = 2;
    private static final int VALUE_COLUMN = 3;
    private static final int ALTERATION_COLUMN = 4;

    // Optional columns
    private static final int GENOTYPE_COLUMN = 5;
    private static final int STATUS_COLUMN = 7;

    private static final String GAIN_LOSS_MISMATCH_WITH_VALUE_ERROR = "gain/loss identifier does not match with value: %s";

    @NotNull
    public static List<CopyNumber> loadCNV(@NotNull final String basePath, @NotNull final String sample) throws IOException, HartwigException {
        final List<CopyNumber> copyNumbers = Lists.newArrayList();
        for (final String line : FreecFileLoader.copyNumberLines(basePath, sample)) {
            copyNumbers.add(FreecCopyNumberFactory.fromCNVLine(line));
        }
        return copyNumbers;
    }

    @NotNull
    public static CopyNumber fromCNVLine(@NotNull final String cnvLine) throws HartwigException {
        final String[] values = cnvLine.split(CNV_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long start = Long.valueOf(values[START_FIELD_COLUMN].trim()) + 1;
        final long end = Long.valueOf(values[END_FIELD_COLUMN].trim());
        final int value = Integer.valueOf(values[VALUE_COLUMN].trim());
        final String alterationString = values[ALTERATION_COLUMN].trim();
        final CopyNumberAlteration alteration = CopyNumberAlteration.fromString(alterationString);

        if (!alteration.equals(CopyNumberAlteration.fromCopyNumber(value))) {
            final String identification = chromosome + "[" + Long.toString(start - 1) + " - " + Long.toString(end) + "]";
            throw new HartwigException(String.format(GAIN_LOSS_MISMATCH_WITH_VALUE_ERROR, identification));
        }

        final ImmutableFreecCopyNumber.Builder builder =
                ImmutableFreecCopyNumber.builder().chromosome(chromosome).start(start).end(end).value(value);

        if (values.length > GENOTYPE_COLUMN) {
            builder.genotype(values[GENOTYPE_COLUMN].trim());
        }

        if (values.length > STATUS_COLUMN) {
            builder.status(FreecCopyNumber.Status.valueOf(values[STATUS_COLUMN].trim().toUpperCase()));
        }

        return builder.build();
    }
}
