package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public enum WriteType
{
    GENERIC,
    TYPE_SPECIFIC,
    SUMMARY;

    public static List<WriteType> DEFAULT_WRITE_TYPES = List.of(GENERIC, TYPE_SPECIFIC);

    public static List<WriteType> fromConfigStr(final String writeTypeStr)
    {
        if(writeTypeStr == null || writeTypeStr.isEmpty())
            return DEFAULT_WRITE_TYPES;

        String[] types = writeTypeStr.split(ITEM_DELIM);
        return Arrays.stream(types).map(x -> WriteType.valueOf(x)).collect(Collectors.toList());
    }
}
