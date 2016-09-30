package com.hartwig.hmftools.ecrfanalyser.reader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CodeListFactory {
    private static final Logger LOGGER = LogManager.getLogger(CodeListFactory.class);

    private static final String CODE_LIST_MARKER = "=";

    private CodeListFactory() {
    }

    @NotNull
    static String fromText(@NotNull String codeListItemString) {
        // KODU: Convert "1=x" to "x".
        int marker = codeListItemString.indexOf(CODE_LIST_MARKER);
        if (marker < 0) {
            LOGGER.warn("Invalid code list item, expected a \"" + CODE_LIST_MARKER + "\" in " + codeListItemString);
            return codeListItemString;
        }

        return codeListItemString.substring(marker + 1).trim();
    }
}
