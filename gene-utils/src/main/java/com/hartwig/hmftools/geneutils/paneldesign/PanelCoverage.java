package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Interface for querying regions already covered by the panel probes.
public interface PanelCoverage
{
    // Checks if a target region is fully covered by probes in the panel.
    boolean isCovered(final ChrBaseRegion target);

    // Gets all regions covered by probes in the panel.
    Stream<ChrBaseRegion> coveredRegions();
}
