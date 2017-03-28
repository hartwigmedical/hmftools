package com.hartwig.hmftools.common.variant;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(GermlineVariantFactory.class);

    private static final String VCF_COLUMN_SEPARATOR = "\t";

    private static final int REF_COLUMN = 3;
    private static final int ALT_COLUMN = 4;
    private static final int FILTER_COLUMN = 6;

    private static final int REF_SAMPLE_COLUMN = 9;
    private static final int TUMOR_SAMPLE_COLUMN = 10;
    private static final int SAMPLE_DATA_GENOTYPE_COLUMN = 0;
    private static final int SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN = 1;

    private GermlineVariantFactory() {
    }

    @Nullable
    public static GermlineVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        if (values.length <= REF_SAMPLE_COLUMN) {
            LOGGER.warn("Not enough columns in vcf line: " + line);
            return null;
        }

        final VariantType type = VariantType.fromRefAlt(values[REF_COLUMN].trim(), values[ALT_COLUMN].trim());
        final String filter = values[FILTER_COLUMN].trim();
        final GermlineSampleData refData = fromSampleData(values[REF_SAMPLE_COLUMN].trim());
        if (refData == null) {
            LOGGER.warn("Could not generate sample data for ref: " + values[REF_SAMPLE_COLUMN].trim());
            return null;
        }

        final GermlineSampleData tumorData =
                values.length > TUMOR_SAMPLE_COLUMN ? fromSampleData(values[TUMOR_SAMPLE_COLUMN].trim()) : null;

        return new GermlineVariant(type, filter, refData, tumorData);
    }

    @Nullable
    @VisibleForTesting
    static GermlineSampleData fromSampleData(@NotNull final String sampleData) {
        final String parts[] = VariantFactoryFunctions.splitSampleDataFields(sampleData);
        if (parts.length < 2) {
            LOGGER.warn("Could not parse germline sample data: " + sampleData);
            return null;
        }

        final AlleleFrequencyData alleleFrequencyData = VariantFactoryFunctions.determineAlleleFrequencies(
                parts[SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN].trim());
        if (alleleFrequencyData == null) {
            LOGGER.warn("Could not parse allele frequencies for germline sample data: " + sampleData);
            return null;
        }

        return new GermlineSampleData(parts[SAMPLE_DATA_GENOTYPE_COLUMN].trim(), alleleFrequencyData.totalReadCount(),
                alleleFrequencyData.alleleReadCount());
    }
}
