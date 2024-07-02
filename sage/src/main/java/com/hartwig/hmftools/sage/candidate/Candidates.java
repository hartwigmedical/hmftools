package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.select.TierSelector;

public class Candidates
{
    private final List<SimpleVariant> mHotspots;
    private final List<BaseRegion> mPanel;
    private final List<BaseRegion> mHighConfidence;
    private Map<SimpleVariant,List<Candidate>> mCandidateMap;
    private final List<Candidate> mCandidateList;

    public Candidates(final List<SimpleVariant> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
    {
        mHotspots = hotspots;
        mPanel = panel;
        mHighConfidence = highConfidence;

        mCandidateMap = null;
        mCandidateList = Lists.newArrayList();
    }

    public void addOfMultipleSamples(final Collection<AltContext> altContexts)
    {
        if(mCandidateMap == null)
            mCandidateMap = Maps.newHashMap();

        final TierSelector tierSelector = new TierSelector(mHotspots, mPanel, mHighConfidence);
        final SimpleVariantComparator variantComparator = new SimpleVariantComparator();

        for(final AltContext altContext : altContexts)
        {
            List<Candidate> candidates = mCandidateMap.get(altContext);

            if(candidates == null)
            {
                candidates = Lists.newArrayList();
                mCandidateMap.put(altContext, candidates);
            }

            Candidate newCandidate = Candidate.fromAltContext(tierSelector.tier(altContext), altContext);

            Candidate matchingCandidate = candidates.stream()
                    .filter(x -> variantComparator.compare(x.variant(), newCandidate.variant()) == 0)
                    .filter(x -> x.readContext().coreStr().equals(newCandidate.readContext().coreStr()))
                    .findFirst().orElse(null);

            if(matchingCandidate != null)
            {
                matchingCandidate.update(altContext);
            }
            else
            {
                candidates.add(newCandidate);
            }
        }
    }

    public void addSingleSample(final Collection<AltContext> altContexts)
    {
        final TierSelector tierSelector = new TierSelector(mHotspots, mPanel, mHighConfidence);

        final SimpleVariantComparator variantComparator = new SimpleVariantComparator();

        for(AltContext altContext : altContexts)
        {
            Candidate candidate = Candidate.fromAltContext(tierSelector.tier(altContext), altContext);

            int index = 0;

            while(index < mCandidateList.size())
            {
                Candidate existingCandidate = mCandidateList.get(index);
                if(variantComparator.compare(existingCandidate.variant(), candidate.variant()) > 0)
                    break;

                ++index;
            }

            mCandidateList.add(index, candidate);
        }
    }

    public List<Candidate> candidates(final List<BasePosition> restrictedPositions)
    {
        if(mCandidateMap != null)
        {
            final SimpleVariantComparator variantComparator = new SimpleVariantComparator();
            mCandidateMap.values().stream().forEach(x -> mCandidateList.addAll(x));
            mCandidateList.sort((o1, o2) -> variantComparator.compare(o1.variant(), o2.variant()));
        }

        if(!restrictedPositions.isEmpty())
        {
            List<Candidate> restrictedList = Lists.newArrayList();

            for(Candidate candidate : mCandidateList)
            {
                if(restrictedPositions.stream().anyMatch(x -> x.matches(candidate.chromosome(), candidate.position())))
                {
                    restrictedList.add(candidate);
                }
            }

            mCandidateList.clear();
            mCandidateList.addAll(restrictedList);
        }

        return mCandidateList;
    }

    private class SimpleVariantComparator implements Comparator<SimpleVariant>
    {
        @Override
        public int compare(final SimpleVariant o1, final SimpleVariant o2)
        {
            int standardCompare = o1.compareTo(o2);

            if(standardCompare != 0)
                return standardCompare;

            int o1Length = Math.max(o1.ref().length(), o1.alt().length());
            int o2Length = Math.max(o2.ref().length(), o2.alt().length());
            int lengthCompare = Integer.compare(o1Length, o2Length);

            if(lengthCompare != 0)
                return lengthCompare;

            int refCompare = o1.ref().compareTo(o2.ref());

            if(refCompare != 0)
                return refCompare;

            return o1.alt().compareTo(o2.alt());
        }
    }

}
