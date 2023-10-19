package com.hartwig.hmftools.sieve.unmap;

public class PrimaryReadInfo
{
    private final String mReadName;
    private final boolean mFirstInPair;

    public PrimaryReadInfo(final String readName, final boolean firstInPair)
    {
        mReadName = readName;
        mFirstInPair = firstInPair;
    }
}
