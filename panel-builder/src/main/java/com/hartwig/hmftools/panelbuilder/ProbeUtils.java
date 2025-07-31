package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.Utils.regionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.Utils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.Utils.regionStartingAt;

import com.hartwig.hmftools.common.region.BaseRegion;

// Miscellaneous probe maths and utilities.
public class ProbeUtils
{
    public static BaseRegion probeRegionStartingAt(int startPosition)
    {
        return regionStartingAt(startPosition, PROBE_LENGTH);
    }

    public static BaseRegion probeRegionCenteredAt(int centrePosition)
    {
        return regionCenteredAt(centrePosition, PROBE_LENGTH);
    }

    public static BaseRegion probeRegionEndingAt(int endPosition)
    {
        return regionEndingAt(endPosition, PROBE_LENGTH);
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
