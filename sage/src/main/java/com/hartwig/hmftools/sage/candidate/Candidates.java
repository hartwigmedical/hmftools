package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.select.TierSelector;

import org.jetbrains.annotations.NotNull;

public class Candidates
{

    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panel;
    private final List<GenomeRegion> highConfidence;
    private final Map<VariantHotspot, Candidate> candidateMap = Maps.newHashMap();
    private List<Candidate> candidateList;

    public Candidates(final List<VariantHotspot> hotspots, final List<GenomeRegion> panel, final List<GenomeRegion> highConfidence)
    {
        this.hotspots = hotspots;
        this.panel = panel;
        this.highConfidence = highConfidence;
    }

    public void add(@NotNull final Collection<AltContext> altContexts)
    {
        if(candidateList != null)
        {
            throw new IllegalStateException("Cannot add more alt contexts");
        }

        final TierSelector tierSelector = new TierSelector(hotspots, panel, highConfidence);
        for(final AltContext altContext : altContexts)
        {
            candidateMap.computeIfAbsent(altContext, x -> new Candidate(tierSelector.tier(altContext), altContext)).update(altContext);
        }
    }

    @NotNull
    public List<Candidate> candidates()
    {
        if(candidateList == null)
        {
            final VariantHotspotComparator variantHotspotComparator = new VariantHotspotComparator();
            candidateList = Lists.newArrayList(candidateMap.values());
            candidateList.sort((o1, o2) -> variantHotspotComparator.compare(o1.variant(), o2.variant()));
        }

        return candidateList;
    }
}
