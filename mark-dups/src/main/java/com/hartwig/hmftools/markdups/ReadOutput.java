package com.hartwig.hmftools.markdups;

import java.util.Arrays;
import java.util.stream.Collectors;

public enum ReadOutput
{
    NONE,
    DUPLICATES,
    MISMATCHES,
    ALL;

    public static String valuesStr()
    {
        return Arrays.stream(ReadOutput.values()).map(x -> x.toString()).collect(Collectors.joining(","));
    }
}
