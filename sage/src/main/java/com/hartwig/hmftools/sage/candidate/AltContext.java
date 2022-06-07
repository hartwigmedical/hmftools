package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS_PERC;

import java.util.Collections;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.ReadContextMatch;

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
    private AltContext mSecondCandidate; // relevant if has a different read context and sufficient support

    public AltContext(final RefContext refContext, final String ref, final String alt)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;

        mReadContextCandidates = Lists.newArrayList();
        mCandidate = null;
        mSecondCandidate = null;
    }

    public AltContext(
            final RefContext refContext, final String ref, final String alt, final ReadContextCandidate candidate,
            int rawSupportAlt, int rawBaseQualAlt)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;

        mCandidate = candidate;
        mRawSupportAlt = rawSupportAlt;
        mRawBaseQualityAlt = rawBaseQualAlt;
        mReadContextCandidates = null;
    }

    public void incrementAltRead(int baseQuality)
    {
        mRawSupportAlt++;
        mRawBaseQualityAlt += baseQuality;
    }

    public void addReadContext(int numberOfEvents, final ReadContext newReadContext)
    {
        int partialMatch = 0;
        int coreMatch = 0;
        ReadContextCandidate fullMatchCandidate = null;

        for(ReadContextCandidate candidate : mReadContextCandidates)
        {
            final ReadContextMatch match = candidate.readContext().indexedBases().matchAtPosition(newReadContext.indexedBases());

            switch(match)
            {
                case FULL:
                    candidate.incrementFull(1, numberOfEvents);
                    fullMatchCandidate = candidate;
                    break;
                case PARTIAL:
                    candidate.PartialMatch++;
                    partialMatch++;
                    break;
                case CORE:
                    candidate.CoreMatch++;
                    coreMatch++;
                    break;
            }
        }

        if(fullMatchCandidate == null)
        {
            final ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, newReadContext);
            candidate.CoreMatch += coreMatch;
            candidate.PartialMatch += partialMatch;
            mReadContextCandidates.add(candidate);
        }
        else if(newReadContext.maxFlankLength() > fullMatchCandidate.maxFlankLength())
        {
            mReadContextCandidates.remove(fullMatchCandidate);
            final ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, newReadContext);
            candidate.CoreMatch += fullMatchCandidate.CoreMatch;
            candidate.PartialMatch += fullMatchCandidate.PartialMatch;
            candidate.incrementFull(fullMatchCandidate.FullMatch, fullMatchCandidate.mMinNumberOfEvents);
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

    public boolean hasValidCandidate() { return mCandidate != null; }
    public boolean hasSecondCandidate() { return mSecondCandidate != null; }
    public AltContext secondCandidate() { return mSecondCandidate; }

    public void selectCandidates()
    {
        if(mReadContextCandidates.isEmpty())
            return;

        mReadContextCandidates.removeIf(x -> x.readContext().hasIncompleteFlanks());

        // sort by full, the partial then core read counts
        Collections.sort(mReadContextCandidates);

        if(mReadContextCandidates.isEmpty())
            return;

        mCandidate = mReadContextCandidates.get(0);

        if(mReadContextCandidates.size() > 1)
        {
            final String topCore = mReadContextCandidates.get(0).readContext().coreString();
            double topCandidateRcThreshold = mReadContextCandidates.get(0).fullMatch() * MIN_SECOND_CANDIDATE_FULL_READS_PERC;

            // add a second if its core is different and it has sufficient support
            for(int i = 0; i < mReadContextCandidates.size(); ++i)
            {
                ReadContextCandidate candidate = mReadContextCandidates.get(i);

                if(candidate.fullMatch() < MIN_SECOND_CANDIDATE_FULL_READS || candidate.fullMatch() < topCandidateRcThreshold)
                    break;

                String coreStr = candidate.readContext().coreString();
                if(coreStr.contains(topCore) || topCore.contains(coreStr))
                    continue;

                mSecondCandidate = new AltContext(RefContext, Ref, Alt, candidate, mRawSupportAlt, mRawBaseQualityAlt);
                break;
            }
        }

        mReadContextCandidates.clear();
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
        return String.format("var(%s:%d %s->%s) readCandidates(%d)",
                chromosome(), position(), Ref, Alt, mReadContextCandidates != null ? mReadContextCandidates.size() : 0);
    }

    protected class ReadContextCandidate implements Comparable<ReadContextCandidate>
    {
        private final ReadContext mReadContext;
        public int FullMatch;
        public int PartialMatch;
        public int CoreMatch;
        public int mMinNumberOfEvents;

        ReadContextCandidate(int numberOfEvents, final ReadContext readContext)
        {
            mReadContext = readContext;
            mMinNumberOfEvents = numberOfEvents;
        }

        public void incrementFull(int count, int numberOfEvents)
        {
            FullMatch += count;
            mMinNumberOfEvents = min(mMinNumberOfEvents, numberOfEvents);
        }

        public int count() { return FullMatch + PartialMatch; }

        public int minNumberOfEvents() { return mMinNumberOfEvents; }

        public int maxFlankLength()
        {
            return mReadContext.maxFlankLength();
        }

        public ReadContext readContext() { return mReadContext; }

        public int fullMatch() { return FullMatch; }

        @Override
        public int compareTo(@NotNull final ReadContextCandidate o)
        {
            int fullCompare = -Integer.compare(FullMatch, o.FullMatch);

            if(fullCompare != 0)
                return fullCompare;

            int partialCompare = -Integer.compare(PartialMatch, o.PartialMatch);

            if(partialCompare != 0)
                return partialCompare;

            return -Integer.compare(CoreMatch, o.CoreMatch);
        }

        public String toString()
        {
            return String.format("matches(f=%d p=%d c=%d)", FullMatch, PartialMatch, CoreMatch);
        }
    }
}
