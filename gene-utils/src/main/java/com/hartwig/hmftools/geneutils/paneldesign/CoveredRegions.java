package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// TODO
// TODO: should this be more like a probe cache to collect all probes?
public class CoveredRegions
{
    private final ArrayList<ChrBaseRegion> mRegions;

    public CoveredRegions()
    {
        mRegions = new ArrayList<>();
    }

    public boolean isCovered(final ChrBaseRegion target)
    {
        return mRegions.stream().anyMatch(target::overlaps);
    }

    public void addRegion(final ChrBaseRegion region)
    {
        mRegions.add(region);
    }

    public void addFromProbe(final EvaluatedProbe probe)
    {
        addRegion(probe.candidate().probeRegion());
    }

    public void addFromProbes(final List<EvaluatedProbe> probes)
    {
        probes.forEach(this::addFromProbe);
    }

    public void clear()
    {
        mRegions.clear();
    }
}
