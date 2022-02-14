package com.hartwig.hmftools.sage.common;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.jetbrains.annotations.Nullable;

public class SageVariant
{
    private final Candidate mCandidate;
    private final Set<String> mFilters;
    private final List<ReadContextCounter> mNormalReadCounters;
    private final List<ReadContextCounter> mTumorReadCounters;

    private int mMixedImpact;

    public SageVariant(
            final Candidate candidate,  final List<ReadContextCounter> normalCounters, final List<ReadContextCounter> tumorReadCounters)
    {
        mCandidate = candidate;
        mNormalReadCounters = normalCounters;
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
                final int[] counts = mTumorReadCounters.get(0).lpsCounts().get(i);
                return counts[0] + counts[1];
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
    public List<int[]> localPhaseSetCounts()
    {
        if(mTumorReadCounters.isEmpty())
            return null;

        return mTumorReadCounters.get(0).lpsCounts();
    }

    public int mixedGermlineImpact() { return mMixedImpact; }
    public void mixedGermlineImpact(final int mixedImpact) { mMixedImpact = mixedImpact; }

    public boolean isPassing() { return mFilters.isEmpty(); }

    public boolean isTumorEmpty() { return mTumorReadCounters.isEmpty(); }
    public boolean isNormalEmpty() { return mNormalReadCounters.isEmpty(); }

    public VariantHotspot variant() { return mCandidate.variant(); }

    public VariantTier tier() { return mCandidate.tier(); }

    public Set<String> filters() { return mFilters; }

    public ReadContext readContext() { return mTumorReadCounters.get(0).readContext(); }

    public List<ReadContextCounter> normalReadCounters() { return mNormalReadCounters; }
    public List<ReadContextCounter> tumorReadCounters() { return mTumorReadCounters; }

    public int totalQuality()
    {
        return mTumorReadCounters.stream().mapToInt(ReadContextCounter::tumorQuality).sum();
    }

    public boolean isIndel()
    {
        return variant().ref().length() != variant().alt().length();
    }

    public boolean isMnv()
    {
        return variant().ref().length() > 1 && variant().ref().length() == variant().alt().length();
    }

    public boolean isSnv()
    {
        return variant().ref().length() == 1 && variant().alt().length() == 1;
    }

    public String toString()
    {
        return String.format("%s filters(%s)",
                mCandidate.toString(), mFilters.toString());
    }
}
