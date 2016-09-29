package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CodeListFactory {
    private static final Logger LOGGER = LogManager.getLogger(CodeListFactory.class);

    private static final String CODE_LIST_MARKER = "=";

    private CodeListFactory() {
    }

    @NotNull
    static Map<Integer, String> extractValuesFromStrings(@NotNull List<String> codeListItems) {
        Map<Integer, String> values = Maps.newHashMap();

        for (String codeListString : codeListItems) {
            int marker = codeListString.indexOf(CODE_LIST_MARKER);
            if (marker < 0) {
                LOGGER.warn("Cannot parse code list string: " + codeListString);
            } else {
                String index = codeListString.substring(0, marker).trim();
                String value = codeListString.substring(marker + 1).trim();

                try {
                    values.put(Integer.valueOf(index), value);
                } catch (NumberFormatException exception) {
                    LOGGER.warn("Could not convert index to Integer: " + index);
                }
            }
        }
        return values;
    }
}
