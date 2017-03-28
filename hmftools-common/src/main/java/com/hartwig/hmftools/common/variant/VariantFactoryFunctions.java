package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class VariantFactoryFunctions {

    private static final String SAMPLE_FIELD_SEPARATOR = ":";
    private static final String ALLELE_FREQUENCY_FIELD_SEPARATOR = ",";

    private VariantFactoryFunctions() {
    }

    @NotNull
    static String[] splitSampleDataFields(@NotNull final String sampleData) {
        return sampleData.split(SAMPLE_FIELD_SEPARATOR);
    }

    @Nullable
    static ReadCount analyzeAlleleFrequencies(@NotNull final String alleleFrequencyField) {
        final String[] afFields = alleleFrequencyField.split(ALLELE_FREQUENCY_FIELD_SEPARATOR);

        if (afFields.length < 2) {
            return null;
        }

        final int readCount = Integer.valueOf(afFields[1]);
        int totalReadCount = 0;
        for (final String afField : afFields) {
            totalReadCount += Integer.valueOf(afField);
        }

        return new ReadCount(readCount, totalReadCount);
    }
}
