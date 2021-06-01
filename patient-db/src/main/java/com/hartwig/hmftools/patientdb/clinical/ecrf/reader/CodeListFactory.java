package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import org.jetbrains.annotations.NotNull;

final class CodeListFactory {

    private static final String CODE_LIST_MARKER = "=";

    private CodeListFactory() {
    }

    @NotNull
    static String fromText(@NotNull String codeListItemString) {
        // Convert "1=x" to "x".
        int marker = codeListItemString.indexOf(CODE_LIST_MARKER);
        if (marker < 0) {
            return codeListItemString;
        }

        return codeListItemString.substring(marker + 1).trim();
    }
}
