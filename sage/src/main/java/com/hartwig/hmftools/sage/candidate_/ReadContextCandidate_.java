package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.min;

import org.jetbrains.annotations.NotNull;

public class ReadContextCandidate_ implements Comparable<ReadContextCandidate_>
{
    public final ReadContext_ ReadContext;
    public int FullMatch;
    public int PartialMatch;
    public int CoreMatch;
    public int MinNumberOfEvents;

    public ReadContextCandidate_(int numberOfEvents, final ReadContext_ readContext)
    {
        ReadContext = readContext;
        MinNumberOfEvents = numberOfEvents;

        FullMatch = 0;
        PartialMatch = 0;
        CoreMatch = 0;
    }

    public void incrementFull(int count, int numberOfEvents)
    {
        FullMatch += count;
        MinNumberOfEvents = min(MinNumberOfEvents, numberOfEvents);
    }

    /**
     * Full + partial matches.
     */
    public int count() { return FullMatch + PartialMatch; }

    /**
     * Sorts via DESC FullMatch, DESC PartialMatch, DESC CoreMatch.
     */
    @Override
    public int compareTo(@NotNull final ReadContextCandidate_ o)
    {
        int fullCompare = -Integer.compare(FullMatch, o.FullMatch);

        if(fullCompare != 0)
            return fullCompare;

        int partialCompare = -Integer.compare(PartialMatch, o.PartialMatch);

        if(partialCompare != 0)
            return partialCompare;

        return -Integer.compare(CoreMatch, o.CoreMatch);
    }

    @Override
    public String toString()
    {
        return String.format("matches(f=%d p=%d c=%d)", FullMatch, PartialMatch, CoreMatch);
    }
}
