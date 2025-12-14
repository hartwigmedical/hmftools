package com.hartwig.hmftools.datamodel.finding.clinicaltranscript;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

final class CsvFileReader {

    private CsvFileReader() {
    }

    @NotNull
    public static Map<String, Integer> getHeadersToDelimiter(@NotNull String fieldsHeader, @NotNull String delimiter) {
        String[] items = fieldsHeader.split(delimiter, -1);
        Map<String, Integer> fieldsIndexMap = Maps.newLinkedHashMap();

        for (int i = 0; i < items.length; ++i) {
            fieldsIndexMap.put(items[i], i);
        }

        return fieldsIndexMap;
    }
}
