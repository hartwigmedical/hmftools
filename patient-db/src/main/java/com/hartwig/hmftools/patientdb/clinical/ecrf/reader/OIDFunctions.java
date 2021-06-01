package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

final class OIDFunctions {

    private static final String OID_PREFIX = "FLD";

    @VisibleForTesting
    static final String OID_SEPARATOR = ".";
    @VisibleForTesting
    static final String NO_CATEGORY = "GENERAL";

    private OIDFunctions() {
    }

    @NotNull
    static String toEcrfFieldName(@NotNull final String OID) {
        String[] fields = OID.split("\\" + OID_SEPARATOR);
        String name = fields[fields.length - 1];
        if (fields.length == 2) {
            return NO_CATEGORY + OID_SEPARATOR + name;
        } else {
            assert fields.length == 3;
            return fields[1] + OID_SEPARATOR + name;
        }
    }

    @NotNull
    static String toOID(@NotNull String ecrfFieldName) {
        String[] fields = ecrfFieldName.split("\\" + OID_SEPARATOR);
        String category = fields[0];
        String name = fields[1];
        if (category.equals(NO_CATEGORY)) {
            return OID_PREFIX + OID_SEPARATOR + name;
        }

        return OID_PREFIX + OID_SEPARATOR + category + OID_SEPARATOR + name;
    }
}
