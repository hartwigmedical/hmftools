package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.List;
import java.util.Map;
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

    private static final int DEBUG_POSITION = -1;
    // private static final int DEBUG_POSITION = 149801;

    public List<Candidate> candidates()
    {
        if(mCandidateList == null)
        {
            final VariantHotspotComparator variantHotspotComparator = new VariantHotspotComparator();
            mCandidateList = Lists.newArrayList(mCandidateMap.values());
            mCandidateList.sort((o1, o2) -> variantHotspotComparator.compare(o1.variant(), o2.variant()));
        }

        if(DEBUG_POSITION > 0)
            mCandidateList = mCandidateList.stream().filter(x -> x.position() == DEBUG_POSITION).collect(Collectors.toList());

        return mCandidateList;
    }
}
