package com.hartwig.hmftools.sage.candidate;

import java.util.Collections;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextMatch;

import org.jetbrains.annotations.NotNull;

public class AltContext implements VariantHotspot
{
    public final String Ref;
    public final String Alt;
    public final RefContext RefContext;
    
    private final List<ReadContextCandidate> mReadContextCandidates;

    private int mRawSupportAlt;
    private int mRawBaseQualityAlt;
    private ReadContextCandidate mCandidate;

    public AltContext(final RefContext refContext, final String ref, final String alt)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;

        mReadContextCandidates = Lists.newArrayList();
    }

    public void incrementAltRead(int baseQuality)
    {
        mRawSupportAlt++;
        mRawBaseQualityAlt += baseQuality;
    }

    public void addReadContext(int numberOfEvents, @NotNull final ReadContext newReadContext)
    {
        if(mCandidate != null)
            throw new IllegalStateException();

        int partialMatch = 0;
        int coreMatch = 0;
        ReadContextCandidate fullMatchCandidate = null;

        for(ReadContextCandidate candidate : mReadContextCandidates)
        {
            final ReadContextMatch match = candidate.readContext().matchAtPosition(newReadContext);
            switch(match)
            {
                case FULL:
                    candidate.incrementFull(1, numberOfEvents);
                    fullMatchCandidate = candidate;
                    break;
                case PARTIAL:
                    candidate.incrementPartial(1);
                    partialMatch++;
                    break;
                case CORE:
                    candidate.incrementCore(1);
                    coreMatch++;
                    break;
            }
        }

        if(fullMatchCandidate == null)
        {
            final ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, newReadContext);
            candidate.incrementCore(coreMatch);
            candidate.incrementPartial(partialMatch);
            mReadContextCandidates.add(candidate);
        }
        else if(newReadContext.maxFlankLength() > fullMatchCandidate.maxFlankLength())
        {
            mReadContextCandidates.remove(fullMatchCandidate);
            final ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, newReadContext);
            candidate.incrementCore(fullMatchCandidate.mCoreMatch);
            candidate.incrementPartial(fullMatchCandidate.mPartialMatch);
            candidate.incrementFull(fullMatchCandidate.mFullMatch, fullMatchCandidate.mMinNumberOfEvents);
            mReadContextCandidates.add(candidate);
        }
    }

    public int readContextSupport()
    {
        return mCandidate.count();
    }

    public int minNumberOfEvents()
    {
        return mCandidate.minNumberOfEvents();
    }

    @VisibleForTesting
    List<ReadContextCandidate> interimReadContexts()
    {
        return mReadContextCandidates;
    }

    public boolean finaliseAndValidate()
    {
        if(mReadContextCandidates.isEmpty())
            return false;

        mReadContextCandidates.removeIf(x -> x.readContext().hasIncompleteFlanks());

        Collections.sort(mReadContextCandidates);

        if(!mReadContextCandidates.isEmpty())
        {
            mCandidate = mReadContextCandidates.get(0);
        }

        mReadContextCandidates.clear();

        return mCandidate != null;
    }

    public ReadContext readContext() { return mCandidate.readContext(); }

    @Override
    public String ref() { return Ref; }

    @Override
    public String alt() { return Alt; }

    @Override
    public String chromosome() { return RefContext.chromosome(); }

    @Override
    public int position() { return RefContext.position(); }

    public int rawAltSupport() { return mRawSupportAlt; }

    public int rawDepth() { return RefContext.rawDepth(); }

    public int rawAltBaseQuality() { return mRawBaseQualityAlt; }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof VariantHotspot && equalTo((VariantHotspot) another);
    }

    private boolean equalTo(final VariantHotspot another)
    {
        return ref().equals(another.ref()) && alt().equals(another.alt()) && chromosome().equals(another.chromosome())
                && position() == another.position();
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + Ref.hashCode();
        h += (h << 5) + Alt.hashCode();
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Ints.hashCode(position());
        return h;
    }

    public String toString()
    {
        return String.format("var(%s:%d %s->%s) readCandidates(%d)", chromosome(), position(), Ref, Alt, mReadContextCandidates.size());
    }

    static class ReadContextCandidate implements Comparable<ReadContextCandidate>
    {
        private final ReadContext mReadContext;
        private int mFullMatch;
        private int mPartialMatch;
        private int mCoreMatch;
        private int mMinNumberOfEvents;

        ReadContextCandidate(int numberOfEvents, final ReadContext readContext)
        {
            mReadContext = readContext.cloneAndExpand();
            mMinNumberOfEvents = numberOfEvents;
        }

        public void incrementFull(int count, int numberOfEvents)
        {
            mFullMatch += count;
            mMinNumberOfEvents = Math.min(mMinNumberOfEvents, numberOfEvents);
        }

        public void incrementPartial(int count)
        {
            mPartialMatch += count;
        }

        public void incrementCore(int count)
        {
            mCoreMatch += count;
        }

        public int count()
        {
            return mFullMatch + mPartialMatch;
        }

        public int minNumberOfEvents() { return mMinNumberOfEvents; }

        public int maxFlankLength()
        {
            return mReadContext.maxFlankLength();
        }

        public ReadContext readContext() { return mReadContext; }

        public int fullMatch()
        {
            return mFullMatch;
        }

        @Override
        public int compareTo(@NotNull final ReadContextCandidate o)
        {
            int fullCompare = -Integer.compare(mFullMatch, o.mFullMatch);

            if(fullCompare != 0)
                return fullCompare;

            int partialCompare = -Integer.compare(mPartialMatch, o.mPartialMatch);

            if(partialCompare != 0)
                return partialCompare;

            return -Integer.compare(mCoreMatch, o.mCoreMatch);
        }

        public String toString()
        {
            return String.format("matches(f=%d p=%d c=%d)", mFullMatch, mPartialMatch, mCoreMatch);
        }
    }
}
