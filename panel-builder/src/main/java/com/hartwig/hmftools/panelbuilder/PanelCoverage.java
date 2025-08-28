package com.hartwig.hmftools.panelbuilder;

import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Interface for querying regions already covered by the panel probes.
public interface PanelCoverage
{
    // Checks if a region is fully covered by probes in the panel.
    boolean isCovered(final ChrBaseRegion region);

    // Gets all regions covered by probes in the panel.
    Stream<ChrBaseRegion> coveredRegions();
}
