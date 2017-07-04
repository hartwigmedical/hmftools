package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.VariantFactory.VCF_COLUMN_SEPARATOR;
import static com.hartwig.hmftools.common.variant.VariantFactory.sampleFromHeaderLine;
import static com.hartwig.hmftools.common.variant.VariantFactory.withLine;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(GermlineVariantFactory.class);

    private static final int REF_SAMPLE_COLUMN = 9;
    private static final int TUMOR_SAMPLE_COLUMN = 10;
    private static final int SAMPLE_DATA_GENOTYPE_COLUMN = 0;
    private static final int SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN = 1;
    private static final int SAMPLE_DATA_COMBINED_DEPTH_COLUMN = 2;

    private GermlineVariantFactory() {
    }

    @NotNull
    public static String refSampleFromHeaderLine(@NotNull final String headerLine) {
        return sampleFromHeaderLine(headerLine, REF_SAMPLE_COLUMN);
    }

    @NotNull
    public static String tumorSampleFromHeaderLine(@NotNull final String headerLine) {
        return sampleFromHeaderLine(headerLine, TUMOR_SAMPLE_COLUMN);
    }

    @VisibleForTesting
    static String[] sortReferenceAndTumor(@NotNull final String first, @NotNull final String second) {
        int intersectionLength = sampleIntersectionLength(first, second);
        if (intersectionLength != 0 && second.length() > intersectionLength) {
            final String secondRemainder = second.substring(intersectionLength);
            if (isBloodOrReference(secondRemainder)) {
                return new String[]{second, first};
            }
        }

        return new String[]{first, second};
    }

    private static boolean isBloodOrReference(String header) {
        return header.toLowerCase().startsWith("bl") || header.toLowerCase().startsWith("r");
    }

    private static int sampleIntersectionLength(@NotNull final String first, @NotNull final String second) {
        int maxIndex = Math.min(first.length(), second.length());

        for (int i = 0; i < maxIndex; i++) {
            if (first.charAt(i) != second.charAt(i)) {
                return i;
            }
        }

        return maxIndex;
    }

    @Nullable
    public static GermlineVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        if (values.length <= REF_SAMPLE_COLUMN) {
            LOGGER.warn("Not enough columns in vcf line: " + line);
            return null;
        }

        final GermlineSampleData refData = fromSampleData(values[REF_SAMPLE_COLUMN].trim());
        if (refData == null) {
            LOGGER.warn("Could not generate sample data for ref: " + values[REF_SAMPLE_COLUMN].trim());
            return null;
        }

        final GermlineSampleData tumorData =
                values.length > TUMOR_SAMPLE_COLUMN ? fromSampleData(values[TUMOR_SAMPLE_COLUMN].trim()) : null;

        final GermlineVariant.Builder builder = ImmutableGermlineVariant.builder().refData(refData).tumorData(tumorData);
        return withLine(builder, values).build();
    }

    @Nullable
    @VisibleForTesting
    static GermlineSampleData fromSampleData(@NotNull final String sampleData) {
        final String parts[] = VariantFactoryFunctions.splitSampleDataFields(sampleData);
        if (parts.length < 3) {
            LOGGER.warn("Could not parse germline sample data: " + sampleData);
            return null;
        }

        final AlleleFrequencyData alleleFrequencyData = VariantFactoryFunctions.determineAlleleFrequencies(
                parts[SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN].trim());
        if (alleleFrequencyData == null) {
            LOGGER.warn("Could not parse allele frequencies for germline sample data: " + sampleData);
            return null;
        }

        return ImmutableGermlineSampleData.of(parts[SAMPLE_DATA_GENOTYPE_COLUMN].trim(),
                alleleFrequencyData.totalReadCount(),
                alleleFrequencyData.alleleReadCount(),
                Integer.valueOf(parts[SAMPLE_DATA_COMBINED_DEPTH_COLUMN]));
    }
}
