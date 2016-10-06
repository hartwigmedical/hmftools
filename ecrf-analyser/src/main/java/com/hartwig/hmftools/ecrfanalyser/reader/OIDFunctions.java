package com.hartwig.hmftools.ecrfanalyser.reader;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

final class OIDFunctions {

    private static final String OID_PREFIX = "FLD";
    private static final String OID_SEPARATOR = ".";

    @VisibleForTesting
    static final String NO_CATEGORY = "GENERAL";

    private OIDFunctions() {
    }

    @NotNull
    static String category(@NotNull String OID) {
        String[] fields = OID.split("\\" + OID_SEPARATOR);
        if (fields.length == 2) {
            return NO_CATEGORY;
        }

        assert fields.length == 3;
        return fields[1];
    }

    @NotNull
    static String fieldName(@NotNull String OID) {
        String[] fields = OID.split("\\" + OID_SEPARATOR);
        return fields[fields.length - 1];
    }

    @NotNull
    static String toOID(@NotNull String category, @NotNull String fieldName) {
        if (category.equals(NO_CATEGORY)) {
            return OID_PREFIX + OID_SEPARATOR + fieldName;
        }
        return OID_PREFIX + OID_SEPARATOR + category + OID_SEPARATOR + fieldName;
    }
}
