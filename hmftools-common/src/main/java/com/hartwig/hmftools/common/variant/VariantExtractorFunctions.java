package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

final class VariantExtractorFunctions {

    private static final int ALT_INDEX = 4;
    private static final int REF_INDEX = 3;
    private static final String MULTIPLE_ALTS_IDENTIFIER = ",";

    private VariantExtractorFunctions() {
    }

    @NotNull
    static VariantType extractVCFType(@NotNull final String[] values) {
        final String refValue = values[REF_INDEX];
        final String altValue = values[ALT_INDEX];

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
