package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

final class VariantExtractorFunctions {

    private static final String MULTIPLE_ALTS_IDENTIFIER = ",";

    private VariantExtractorFunctions() {
    }

    @NotNull
    static VariantType extractVCFType(@NotNull final String refValue, @NotNull final String altValue) {
        final String[] allAlts = altValue.split(MULTIPLE_ALTS_IDENTIFIER);

        VariantType type = VariantType.SNP;

        for (String alt : allAlts) {
            if (refValue.length() != alt.length()) {
                type = VariantType.INDEL;
            }
        }
        return type;
    }
}
