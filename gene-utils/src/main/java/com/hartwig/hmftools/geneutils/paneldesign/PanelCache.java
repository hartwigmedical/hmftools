package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PanelCache
{
    private final Map<String,List<PanelRegion>> mChrRegionsMap;

    public PanelCache()
    {
        mChrRegionsMap = Maps.newHashMap();
    }

    public Map<String,List<PanelRegion>> chrRegionsMap() { return mChrRegionsMap; }

    public List<PanelRegion> chromosomeRegions(final String chromosome)
    {
        return mChrRegionsMap.containsKey(chromosome) ? mChrRegionsMap.get(chromosome) : Collections.emptyList();
    }

    public int regionCount() { return mChrRegionsMap.values().stream().mapToInt(x -> x.size()).sum(); }

    public void addRegion(final PanelRegion panelRegion) { addRegion(panelRegion, true); }

    public void addRegion(final PanelRegion panelRegion, boolean checkOverlaps)
    {
        // don't merge or order regions at this point

        List<PanelRegion> regions = mChrRegionsMap.get(panelRegion.Chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrRegionsMap.put(panelRegion.Chromosome, regions);
        }

        if(checkOverlaps)
        {
            PanelRegion existing = regions.stream().filter(x -> x.overlaps(panelRegion)).findFirst().orElse(null);

            if(existing != null)
            {
                GU_LOGGER.warn("panel region({}) overlaps with existing({})", panelRegion, existing);
            }
        }

        regions.add(panelRegion);
    }

    public boolean overlapsExisting(final String chromosome, final int positionStart, final int positionEnd)
    {
        List<PanelRegion> regions = mChrRegionsMap.get(chromosome);
        return regions != null && regions.stream().anyMatch(x -> x.overlaps(chromosome, positionStart, positionEnd));
    }

    public boolean overlapsExisting(final ChrBaseRegion region)
    {
        List<PanelRegion> regions = mChrRegionsMap.get(region.Chromosome);
        return regions != null && regions.stream().anyMatch(x -> x.overlaps(region));
    }
}
