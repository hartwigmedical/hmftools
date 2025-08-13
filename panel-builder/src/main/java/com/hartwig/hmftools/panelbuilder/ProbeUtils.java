package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Miscellaneous probe maths and utilities.
public class ProbeUtils
{
    public static BaseRegion probeRegionStartingAt(int startPosition)
    {
        return regionStartingAt(startPosition, PROBE_LENGTH);
    }

    public static ChrBaseRegion probeRegionStartingAt(final String chromosome, int startPosition)
    {
        return ChrBaseRegion.from(chromosome, probeRegionStartingAt(startPosition));
    }

    public static BaseRegion probeRegionCenteredAt(int centrePosition)
    {
        return regionCenteredAt(centrePosition, PROBE_LENGTH);
    }

    public static ChrBaseRegion probeRegionCenteredAt(final BasePosition centrePosition)
    {
        return probeRegionCenteredAt(centrePosition.Chromosome, centrePosition.Position);
    }

    public static ChrBaseRegion probeRegionCenteredAt(final String chromosome, int centrePosition)
    {
        return ChrBaseRegion.from(chromosome, probeRegionCenteredAt(centrePosition));
    }

    public static BaseRegion probeRegionEndingAt(int endPosition)
    {
        return regionEndingAt(endPosition, PROBE_LENGTH);
    }

    public static ChrBaseRegion probeRegionEndingAt(final String chromosome, int endPosition)
    {
        return ChrBaseRegion.from(chromosome, probeRegionEndingAt(endPosition));
    }

    // Calculates the minimum probe starting position such that the specified position is contained within the probe.
    public static int minProbeStartContaining(int targetPosition)
    {
        return probeRegionEndingAt(targetPosition).start();
    }

    // Calculates the maximum probe ending position such that the specified position is contained within the probe.
    public static int maxProbeEndContaining(int targetPosition)
    {
        return probeRegionStartingAt(targetPosition).end();
    }

    // Calculates the minimum probe starting position such that the specified region overlaps the probe.
    public static int minProbeStartOverlapping(final BaseRegion region)
    {
        return minProbeStartContaining(region.start());
    }

    // Calculates the maximum probe ending position such that the specified region overlaps the probe.
    public static int maxProbeEndOverlapping(final BaseRegion region)
    {
        return maxProbeEndContaining(region.end());
    }

    // Calculates the minimum probe starting position such that the specified region is overlapping or directly adjacent to the probe.
    public static int minProbeStartWithoutGap(final BaseRegion region)
    {
        return minProbeStartOverlapping(region) - 1;
    }

    // Calculates the maximum probe ending position such that the specified region is overlapping or directly adjacent to the probe.
    public static int maxProbeEndWithoutGap(final BaseRegion region)
    {
        return maxProbeEndOverlapping(region) + 1;
    }
}
