package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class PanelCache
{
    private final Map<String, List<ProbeCandidate>> mChrRegionsMap;

    public PanelCache()
    {
        mChrRegionsMap = Maps.newHashMap();
    }

    public Map<String, List<ProbeCandidate>> chrRegionsMap()
    {
        return mChrRegionsMap;
    }

    public List<ProbeCandidate> chromosomeRegions(final String chromosome)
    {
        return mChrRegionsMap.getOrDefault(chromosome, Collections.emptyList());
    }

    public void addRegion(final ProbeCandidate panelRegion)
    {
        addRegion(panelRegion, true);
    }

    public void addRegion(final ProbeCandidate panelRegion, boolean checkOverlaps)
    {
        // don't merge or order regions at this point

        List<ProbeCandidate> regions = mChrRegionsMap.get(panelRegion.Chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrRegionsMap.put(panelRegion.Chromosome, regions);
        }

        if(checkOverlaps)
        {
            regions.stream()
                    .filter(x -> x.overlaps(panelRegion))
                    .findFirst()
                    .ifPresent(existing -> GU_LOGGER.warn("panel region({}) overlaps with existing({})", panelRegion, existing));
        }

        regions.add(panelRegion);
    }

    public boolean overlapsExisting(final String chromosome, final int positionStart, final int positionEnd)
    {
        List<ProbeCandidate> regions = mChrRegionsMap.get(chromosome);
        return regions != null && regions.stream().anyMatch(x -> x.overlaps(chromosome, positionStart, positionEnd));
    }
}
