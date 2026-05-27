package com.hartwig.hmftools.linx.cohort;

import static java.lang.String.format;

public class PonMatch
{
    public final PonMatchType Type;
    public final int Count;

    public static final PonMatch NONE = new PonMatch(PonMatchType.NONE, 0);

    public PonMatch(final PonMatchType type, final int count)
    {
        Type = type;
        Count = count;
    }

    public String toString() { return format("%s=%d", Type, Count); }
}
