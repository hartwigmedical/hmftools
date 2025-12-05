package com.hartwig.hmftools.panelbuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Implements PanelCoverage for localised checks of probe coverage.
public class SimpleCoverageBuffer implements PanelCoverage
{
    private final List<ChrBaseRegion> mCoveredRegions = new ArrayList<>();

    public void addCoveredRegion(final ChrBaseRegion region)
    {
        mCoveredRegions.add(region);
    }

    public void addCoveredRegions(final List<ChrBaseRegion> regions)
    {
        mCoveredRegions.addAll(regions);
    }

    public void addCoveredRegions(Stream<ChrBaseRegion> regions)
    {
        regions.forEach(this::addCoveredRegion);
    }

    public void clear()
    {
        mCoveredRegions.clear();
    }

    @Override
    public Stream<ChrBaseRegion> coveredRegions()
    {
        return mCoveredRegions.stream();
    }
}
