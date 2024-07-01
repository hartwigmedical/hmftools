package com.hartwig.hmftools.sage.common;

import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.LONG_GERMLINE_INSERT_LENGTH;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.filter.SoftFilter;

import org.jetbrains.annotations.Nullable;

public class SageVariant
{
    private final Candidate mCandidate;
    private final Set<SoftFilter> mFilters;
    private final List<ReadContextCounter> mReferenceReadCounters;
    private final List<ReadContextCounter> mTumorReadCounters;

    private int mMixedImpact;

    public SageVariant(
            final Candidate candidate,  final List<ReadContextCounter> referenceCounters, final List<ReadContextCounter> tumorReadCounters)
    {
        mCandidate = candidate;
        mReferenceReadCounters = referenceCounters;
        mTumorReadCounters = tumorReadCounters;
        mFilters = Sets.newHashSet();
    }

    public Candidate candidate()
    {
        return mCandidate;
    }

    public String chromosome()
    {
        return variant().chromosome();
    }
    public int position() { return variant().position(); }

    public String ref()
    {
        return variant().ref();
    }
    public String alt()
    {
        return variant().alt();
    }

    public int end()
    {
        return position() + ref().length() - 1;
    }

    @Nullable
    public List<Integer> localPhaseSets()
    {
        if(mTumorReadCounters.isEmpty())
            return null;

        return mTumorReadCounters.get(0).localPhaseSets();
    }

    public boolean hasLocalPhaseSets()
    {
        if(mTumorReadCounters.isEmpty())
            return false;

        return mTumorReadCounters.get(0).localPhaseSets() != null && !mTumorReadCounters.get(0).localPhaseSets().isEmpty();
    }

    public boolean hasMatchingLps(final List<Integer> otherLocalPhaseSets)
    {
        if(otherLocalPhaseSets == null)
            return false;

        if(!hasLocalPhaseSets())
            return false;

        return mTumorReadCounters.get(0).localPhaseSets().stream().anyMatch(x -> otherLocalPhaseSets.contains(x));
    }

    public boolean hasMatchingLps(final Set<Integer> otherLocalPhaseSets)
    {
        if(otherLocalPhaseSets == null)
            return false;

        if(!hasLocalPhaseSets())
            return false;

        return mTumorReadCounters.get(0).localPhaseSets().stream().anyMatch(x -> otherLocalPhaseSets.contains(x));
    }

    public boolean hasMatchingLps(final Integer lps)
    {
        if(!hasLocalPhaseSets())
            return false;

        return mTumorReadCounters.get(0).localPhaseSets().contains(lps);
    }

    public int getLpsReadCount(int lps)
    {
        if(!hasLocalPhaseSets())
            return 0;

        for(int i = 0; i < mTumorReadCounters.get(0).localPhaseSets().size(); ++i)
        {
            if(mTumorReadCounters.get(0).localPhaseSets().get(i) == lps)
            {
                return mTumorReadCounters.get(0).lpsCounts().get(i);
            }
        }

        return 0;
    }

    public void removeLps(int lps)
    {
        if(!hasLocalPhaseSets())
            return;

        for(int i = 0; i < mTumorReadCounters.get(0).localPhaseSets().size(); ++i)
        {
            if(mTumorReadCounters.get(0).localPhaseSets().get(i) == lps)
            {
                mTumorReadCounters.get(0).localPhaseSets().remove(i);
                mTumorReadCounters.get(0).lpsCounts().remove(i);
                return;
            }
        }
    }

    @Nullable
    public List<Integer> localPhaseSetCounts()
    {
        if(mTumorReadCounters.isEmpty())
            return null;

        return mTumorReadCounters.get(0).lpsCounts();
    }

    public int mixedGermlineImpact() { return mMixedImpact; }
    public void mixedGermlineImpact(final int mixedImpact) { mMixedImpact = mixedImpact; }

    public boolean isPassing() { return mFilters.isEmpty(); }

    public boolean hasTumorSamples() { return !mTumorReadCounters.isEmpty(); }
    public boolean hasReferenceSamples() { return !mReferenceReadCounters.isEmpty(); }

    public SimpleVariant variant() { return mCandidate.variant(); }

    public VariantTier tier() { return mCandidate.tier(); }

    public Set<SoftFilter> filters() { return mFilters; }
    public Set<String> filtersStringSet() { return mFilters.stream().map(x -> x.filterName()).collect(Collectors.toSet()); }
    public String filtersStr() { return mFilters.stream().map(x -> x.filterName()).collect(Collectors.joining(",")); }

    public VariantReadContext readContext() { return mTumorReadCounters.get(0).readContext(); }

    public List<ReadContextCounter> referenceReadCounters() { return mReferenceReadCounters; }
    public List<ReadContextCounter> tumorReadCounters() { return mTumorReadCounters; }

    public int totalQuality() { return mTumorReadCounters.stream().mapToInt(x -> (int)round(x.tumorQuality())).sum(); }

    public boolean isIndel() { return variant().ref().length() != variant().alt().length(); }
    public boolean isMnv()
    {
        return variant().ref().length() > 1 && variant().ref().length() == variant().alt().length();
    }
    public boolean isSnv()
    {
        return variant().ref().length() == 1 && variant().alt().length() == 1;
    }
    public boolean isDelete() { return variant().ref().length() > variant().alt().length(); }
    public boolean isInsert() { return variant().ref().length() < variant().alt().length(); }

    public boolean isLongInsert() { return SimpleVariant.isLongInsert(mCandidate.variant()); }

    public String toString()
    {
        return String.format("%s filters(%s)",
                mCandidate.toString(), mFilters.toString());
    }
}
