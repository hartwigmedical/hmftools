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
        FullMatch = 0;
        CoreMatch = 0;
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

    public String chromosome() { return mReadContext.variant().Chromosome; }
    public int position() { return mReadContext.variant().Position; }
    public String ref() { return mReadContext.ref(); }
    public String alt() { return mReadContext.alt(); }

    @Override
    public int compareTo(@NotNull final ReadContextCandidate other)
    {
        int fullCompare = -Integer.compare(FullMatch, other.FullMatch);

        if(fullCompare != 0)
            return fullCompare;

        int coreCompare = -Integer.compare(CoreMatch, other.CoreMatch);

        if(coreCompare != 0)
            return coreCompare;

        return Integer.compare(MinNumberOfEvents, other.MinNumberOfEvents);
    }

    public String toString()
    {
        return String.format("matches(full=%d core=%d)", FullMatch, CoreMatch);
    }
}
