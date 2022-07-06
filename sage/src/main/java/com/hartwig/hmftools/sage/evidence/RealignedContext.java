package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class RealignedContext
{
    public final RealignedType Type;
    public final int MatchLength;
    public final int MatchReadIndex;
    public final int RepeatCount;
    public final int RepeatLength;

    public static final RealignedContext NONE = new RealignedContext(
            RealignedType.NONE, 0, 0, 0, 0);

    public RealignedContext(
            final RealignedType type, final int matchLength, final int matchReadIndex, final int repeatCount, final int repeatLength)
    {
        Type = type;
        MatchLength = matchLength;
        MatchReadIndex = matchReadIndex;
        RepeatCount = repeatCount;
        RepeatLength = repeatLength;
    }

    public String toString() { return format("%s match(len=%d index=%d) repeat(%d len=%d)",
            Type, MatchLength, MatchReadIndex, RepeatCount, RepeatLength); }
}
