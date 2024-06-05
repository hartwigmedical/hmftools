package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import java.util.StringJoiner;

// purely used for statistics
public enum ReadMatchType
{
    ALT_SUPPORT,
    REF_SUPPORT,
    SOFT_CLIP,
    CHIMERIC,
    MAX_COVERAGE,
    MAP_QUAL,
    BASE_QUAL,
    NON_CORE,
    IN_SPLIT,
    UNRELATED;

    public static String countsToString(final long[] counts)
    {
        StringJoiner sj = new StringJoiner(", ");
        for(ReadMatchType type : ReadMatchType.values())
        {
            long count = counts[type.ordinal()];
            if(count > 0)
                sj.add(format("%s=%d", type, count));
        }

        return sj.toString();
    }
}
