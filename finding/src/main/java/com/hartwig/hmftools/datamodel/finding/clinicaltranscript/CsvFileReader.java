package com.hartwig.hmftools.datamodel.finding.clinicaltranscript;

import java.util.LinkedHashMap;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

final class CsvFileReader {

    private CsvFileReader() {
    }

    @NotNull
    public static Map<String, Integer> getHeadersToDelimiter(@NotNull String fieldsHeader, @NotNull String delimiter) {
        String[] items = fieldsHeader.split(delimiter, -1);
        Map<String, Integer> fieldsIndexMap = new LinkedHashMap<>();

        for (int i = 0; i < items.length; ++i) {
            fieldsIndexMap.put(items[i], i);
        }

        return fieldsIndexMap;
    }
}
