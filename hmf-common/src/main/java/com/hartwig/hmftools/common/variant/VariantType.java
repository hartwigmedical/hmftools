package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public enum VariantType {
    MNP,
    SNP,
    INDEL,
    UNDEFINED;

    private static final String MULTIPLE_ALTS_IDENTIFIER = ",";

    @NotNull
    public static VariantType fromRefAlt(@NotNull final String refValue, @NotNull final String altValue) {
        final String[] allAlts = altValue.split(MULTIPLE_ALTS_IDENTIFIER);

        VariantType type = VariantType.SNP;

        for (final String alt : allAlts) {
            if (refValue.length() != alt.length()) {
                type = VariantType.INDEL;
            }
        }
        return type;
    }
}
