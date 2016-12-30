package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(GermlineVariantFactory.class);

    private static final String VCF_COLUMN_SEPARATOR = "\t";

    private static final int REF_COLUMN = 3;
    private static final int ALT_COLUMN = 4;

    private static final int REF_SAMPLE_COLUMN = 9;
    private static final int TUMOR_SAMPLE_COLUMN = 10;

    private GermlineVariantFactory() {
    }

    @Nullable
    public static GermlineVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        if (values.length <= TUMOR_SAMPLE_COLUMN) {
            LOGGER.warn("Invalid vcf line: " + line);
            return null;
        }

        final VariantType type = VariantExtractorFunctions.extractVCFType(values[REF_COLUMN].trim(),
                values[ALT_COLUMN].trim());
        final String refData = values[REF_SAMPLE_COLUMN];
        final String tumData = values[TUMOR_SAMPLE_COLUMN];

        return new GermlineVariant(type, refData, tumData);
    }
}
