package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.min;

import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.NotNull;

public class ReadContextCandidate implements Comparable<ReadContextCandidate>
{
    private final VariantReadContext mReadContext;
    private final ReadContextMatcher mMatcher;

    public int FullMatch;
    public int CoreMatch;
    public int MinNumberOfEvents;
    public int LowQualInCoreCount;

    ReadContextCandidate(int numberOfEvents, final VariantReadContext readContext)
    {
        mReadContext = readContext;
        mMatcher = new ReadContextMatcher(mReadContext, false, false);
        MinNumberOfEvents = numberOfEvents;
        LowQualInCoreCount = 0;
    }

    public void incrementFull(int count, int numberOfEvents)
    {
        FullMatch += count;
        MinNumberOfEvents = min(MinNumberOfEvents, numberOfEvents);
    }

    public VariantReadContext readContext()
    {
        return mReadContext;
    }
    public ReadContextMatcher matcher()
    {
        return mMatcher;
    }

    @Override
    public int compareTo(@NotNull final ReadContextCandidate other)
    {
        int fullCompare = -Integer.compare(FullMatch, other.FullMatch);

        if(fullCompare != 0)
        {
            return fullCompare;
        }

        return -Integer.compare(CoreMatch, other.CoreMatch);
    }

    public String toString()
    {
        return String.format("matches(full=%d core=%d)", FullMatch, CoreMatch);
    }
}
