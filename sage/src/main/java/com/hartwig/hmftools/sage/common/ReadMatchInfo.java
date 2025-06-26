package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

public class ReadMatchInfo
{
    public final ReadContextMatch MatchType;
    public final boolean ExactMatch;

    public static final ReadMatchInfo NO_MATCH = new ReadMatchInfo(ReadContextMatch.NONE, false);

    public ReadMatchInfo(final ReadContextMatch matchType, final boolean exactMatch)
    {
        MatchType = matchType;
        ExactMatch = exactMatch;
    }

    public String toString()
    {
        if(MatchType == ReadContextMatch.NONE)
            return ReadContextMatch.NONE.toString();

        return format("%s %s", MatchType, ExactMatch ? "exact" : "low-qual mismatches");
    }
}
