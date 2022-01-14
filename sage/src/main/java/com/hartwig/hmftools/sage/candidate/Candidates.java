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
    private final Map<VariantHotspot,Candidate> mCandidateMap;
    private List<Candidate> mCandidateList;

    public Candidates(final List<VariantHotspot> hotspots, final List<BaseRegion> panel, final List<BaseRegion> highConfidence)
    {
        mHotspots = hotspots;
        mPanel = panel;
        mHighConfidence = highConfidence;

        mCandidateMap = Maps.newHashMap();
    }

    public void add(final Collection<AltContext> altContexts)
    {
        if(mCandidateList != null)
        {
            throw new IllegalStateException("Cannot add more alt contexts");
        }

        final TierSelector tierSelector = new TierSelector(mHotspots, mPanel, mHighConfidence);

        for(final AltContext altContext : altContexts)
        {
            Candidate candidate = mCandidateMap.get(altContext);

            if(candidate == null)
            {
                candidate = Candidate.fromAltContext(tierSelector.tier(altContext), altContext);
                mCandidateMap.put(altContext, candidate);
            }

            candidate.update(altContext);
        }
    }

    public List<Candidate> candidates(final Set<Integer> restrictedPositions)
    {
        if(mCandidateList == null)
        {
            final VariantHotspotComparator variantHotspotComparator = new VariantHotspotComparator();
            mCandidateList = Lists.newArrayList(mCandidateMap.values());
            mCandidateList.sort((o1, o2) -> variantHotspotComparator.compare(o1.variant(), o2.variant()));
        }

        if(!restrictedPositions.isEmpty())
        {
            mCandidateList = mCandidateList.stream()
                    .filter(x -> restrictedPositions.contains(x.position())).collect(Collectors.toList());
        }

        return mCandidateList;
    }
}
