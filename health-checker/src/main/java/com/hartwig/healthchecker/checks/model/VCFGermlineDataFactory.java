package com.hartwig.healthchecker.checks.model;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VCFGermlineDataFactory {

    private static final Logger LOGGER = LogManager.getLogger(VCFGermlineDataFactory.class);

    private static final String VCF_COLUMN_SEPARATOR = "\t";
    private static final int TUMOR_SAMPLE_COLUMN = 10;
    private static final int REF_SAMPLE_COLUMN = 9;

    private VCFGermlineDataFactory() {
    }

    @Nullable
    public static VCFGermlineData fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        if (values.length <= TUMOR_SAMPLE_COLUMN) {
            LOGGER.warn("Invalid vcf: " + line);
            return null;
        }

        final VCFType type = VCFExtractorFunctions.extractVCFType(values);
        final String refData = values[REF_SAMPLE_COLUMN];
        final String tumData = values[TUMOR_SAMPLE_COLUMN];

        return new VCFGermlineData(type, refData, tumData);
    }
}
