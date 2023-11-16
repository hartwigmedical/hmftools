package com.hartwig.hmftools.sage.sync;

public class FragmentData
{
    public final int FirstPosStart;
    public final int FirstPosEnd;
    public final int SecondPosStart;
    public final int SecondPosEnd;

    public FragmentData(final int firstPosStart, final int firstPosEnd, final int secondPosStart, final int secondPosEnd)
    {
        FirstPosStart = firstPosStart;
        FirstPosEnd = firstPosEnd;
        SecondPosStart = secondPosStart;
        SecondPosEnd = secondPosEnd;
    }
}
