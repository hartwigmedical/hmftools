package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public enum WriteType
{
    SV_DATA,
    SAMPLE_SUMMARY,
    BREAKEND,
    PON;

    public static List<WriteType> DEFAULT_WRITE_TYPES = List.of(SAMPLE_SUMMARY, BREAKEND);

    public static List<WriteType> fromConfigStr(final String writeTypeStr)
    {
        if(writeTypeStr == null || writeTypeStr.isEmpty())
            return DEFAULT_WRITE_TYPES;

        String[] types = writeTypeStr.split(ITEM_DELIM);
        return Arrays.stream(types).map(x -> WriteType.valueOf(x)).collect(Collectors.toList());
    }

}
