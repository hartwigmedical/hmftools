package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

final class CodeListFactory {

    private static final String CODE_LIST_MARKER = "=";
    private CodeListFactory() {
    }

    @NotNull
    static Map<Integer, String> extractValuesFromStrings(@NotNull List<String> codeListItems) {
        Map<Integer, String> values = Maps.newHashMap();

        for (String codeListString : codeListItems) {
            int marker = codeListString.indexOf(CODE_LIST_MARKER);
            String index = codeListString.substring(0, marker).trim();
            String value = codeListString.substring(marker + 1).trim();
            values.put(Integer.valueOf(index), value);
        }
        return values;
    }
}
