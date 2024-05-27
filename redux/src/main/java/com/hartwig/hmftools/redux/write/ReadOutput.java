package com.hartwig.hmftools.redux.write;

import java.util.Arrays;
import java.util.stream.Collectors;

public enum ReadOutput
{
    NONE, // no verbose read data
    DUPLICATES, // any duplicate or primary read
    MISMATCHES, // differences between input BAM read duplicate status (ie from Sambamba) and internally dertermied
    ALL; // all reads written verbosely

    public static String valuesStr()
    {
        return Arrays.stream(ReadOutput.values()).map(x -> x.toString()).collect(Collectors.joining(","));
    }
}
