package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.select.TierSelector;

public class Candidates
{
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanel;
    private final List<BaseRegion> mHighConfidence;
    private Map<VariantHotspot,List<Candidate>> mCandidateMap;
    private final List<Candidate> mCandidateList;

    public Candidates(final List<VariantHotspot> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
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
        final VariantHotspotComparator variantComparator = new VariantHotspotComparator();

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
                    .filter(x -> x.readContext().coreString().equals(newCandidate.readContext().coreString()))
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

        final VariantHotspotComparator variantComparator = new VariantHotspotComparator();

        for(final AltContext altContext : altContexts)
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

    public List<Candidate> candidates(final Set<Integer> restrictedPositions)
    {
        if(mCandidateMap != null)
        {
            final VariantHotspotComparator variantComparator = new VariantHotspotComparator();
            mCandidateMap.values().stream().forEach(x -> mCandidateList.addAll(x));
            mCandidateList.sort((o1, o2) -> variantComparator.compare(o1.variant(), o2.variant()));
        }

        if(!restrictedPositions.isEmpty())
        {
            List<Candidate> restrictedList = mCandidateList.stream()
                    .filter(x -> restrictedPositions.contains(x.position())).collect(Collectors.toList());

            mCandidateList.clear();
            mCandidateList.addAll(restrictedList);
        }

        return mCandidateList;
    }
}
