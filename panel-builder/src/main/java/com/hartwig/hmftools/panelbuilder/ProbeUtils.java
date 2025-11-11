package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.getSubregion;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
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

    public static List<ChrBaseRegion> probeTargetedRegions(final SequenceDefinition definition, final TargetedRange targetedRange)
    {
        // Determine subset of probe regions defined by targetedRange.

        ArrayList<ChrBaseRegion> targetedRegions = new ArrayList<>(2);

        int probeLength = definition.baseLength();
        int targetedStart = targetedRange.startOffset();
        int targetedEnd = targetedRange.endOffset();

        // Compute the intersection with the start region.
        ChrBaseRegion startRegion = definition.startRegion();
        if(startRegion != null)
        {
            int startOffset = 0;
            int endOffset = startRegion.baseLength();
            int intersectionStart = max(startOffset, targetedStart);
            int intersectionEnd = min(endOffset, targetedEnd);
            if(intersectionStart < intersectionEnd)
            {
                Orientation orientation = definition.startOrientation() == null ? Orientation.FORWARD : definition.startOrientation();
                targetedRegions.add(getSubregion(startRegion, orientation, intersectionStart, intersectionEnd));
            }
        }

        // Compute the intersection with the end region.
        ChrBaseRegion endRegion = definition.endRegion();
        if(endRegion != null)
        {
            int startOffset = probeLength - endRegion.baseLength();
            int endOffset = probeLength;
            int intersectionStart = max(startOffset, targetedStart) - startOffset;
            int intersectionEnd = min(endOffset, targetedEnd) - startOffset;
            if(intersectionStart < intersectionEnd)
            {
                Orientation orientation = definition.endOrientation() == null ? Orientation.FORWARD : definition.endOrientation();
                targetedRegions.add(getSubregion(endRegion, orientation, intersectionStart, intersectionEnd));
            }
        }

        return targetedRegions;
    }
}
