package com.hartwig.hmftools.common.ecrf.reader;

import org.jetbrains.annotations.NotNull;

final class CodeListFactory {

    private static final String CODE_LIST_MARKER = "=";

    private CodeListFactory() {
    }

    @NotNull
    static String fromText(@NotNull String codeListItemString) {
        // KODU: Convert "1=x" to "x".
        final int marker = codeListItemString.indexOf(CODE_LIST_MARKER);
        if (marker < 0) {
            return codeListItemString;
        }

        return codeListItemString.substring(marker + 1).trim();
    }
}
