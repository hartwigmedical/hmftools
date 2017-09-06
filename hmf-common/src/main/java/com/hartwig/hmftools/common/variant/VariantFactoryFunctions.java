package com.hartwig.hmftools.common.variant;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;

public final class VariantFactoryFunctions {

    private static final String SAMPLE_FIELD_SEPARATOR = ":";
    private static final String ALLELE_FREQUENCY_FIELD_SEPARATOR = ",";

    private VariantFactoryFunctions() {
    }

    @NotNull
    static String[] splitSampleDataFields(@NotNull final String sampleData) {
        return sampleData.split(SAMPLE_FIELD_SEPARATOR);
    }

    @Nullable
    static AlleleFrequencyData determineAlleleFrequencies(@NotNull final String alleleFrequencyField) {
        final String[] afFields = alleleFrequencyField.split(ALLELE_FREQUENCY_FIELD_SEPARATOR);

        if (afFields.length < 2) {
            return null;
        }

        final int alleleReadCount = Integer.valueOf(afFields[1]);
        int totalReadCount = 0;
        for (final String afField : afFields) {
            totalReadCount += Integer.valueOf(afField);
        }

        return new AlleleFrequencyData(alleleReadCount, totalReadCount);
    }

    @NotNull
    public static AlleleFrequencyData determineAlleleFrequencies(@NotNull final Genotype genotype) {
        Preconditions.checkArgument(genotype.hasAD());

        int[] adFields = genotype.getAD();
        int totalReadCount = 0;
        final int alleleReadCount = adFields[1];
        for (final int afField : adFields) {
            totalReadCount += afField;
        }

        return new AlleleFrequencyData(alleleReadCount, totalReadCount);
    }
}
