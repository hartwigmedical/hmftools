package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

public class RepeatContext
{
    public final byte[] Bases;
    public final int RepeatIndex;
    public final int StartIndex;
    public final int EndIndex;
    public final int Length;
    public final int ForwardCount;
    public final int BackwardCount;

    public RepeatContext(
            final byte[] bases, final int repeatIndex, final int startIndex, final int endIndex, final int length,
            int forwardCount, int backwardCount)
    {
        Bases = bases;
        RepeatIndex = repeatIndex;
        Length = length;
        BackwardCount = backwardCount;
        ForwardCount = forwardCount;
        StartIndex = startIndex;
        EndIndex = endIndex;
    }

    public int count()
    {
        return ForwardCount + BackwardCount;
    }

    public String sequence() { return Length > 0 ? new String(Bases, RepeatIndex, Length) : ""; }

    public String toString()
    {
        if(Length == 0) return "";

        return format("bases(%d len=%d) index(%d start=%d end=%d) repeats(fwd=%d back=%d",
                new String(Bases, RepeatIndex, Length), Length, RepeatIndex, StartIndex, EndIndex, ForwardCount, BackwardCount);
    }
}