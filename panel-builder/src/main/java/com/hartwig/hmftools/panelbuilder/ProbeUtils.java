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

    public static List<TargetRegion> probeTargetRegions(final Probe probe)
    {
        return probeTargetedRegions(probe.definition(), probe.targetedRange()).stream()
                .map(region -> new TargetRegion(region, probe.metadata()))
                .toList();
    }

    // Determines the subset of the probe's reference genome regions covered by targetedRange (offsets into the probe sequence).
    // Handles any number of segments; insert segments occupy sequence space but map to no genome region.
    public static List<ChrBaseRegion> probeTargetedRegions(final SequenceDefinition definition, final TargetedRange targetedRange)
    {
        int targetedStart = targetedRange.startOffset();
        int targetedEnd = targetedRange.endOffset();

        ArrayList<ChrBaseRegion> targetedRegions = new ArrayList<>(definition.segments().size());
        int segmentStart = 0;
        for(SequenceSegment segment : definition.segments())
        {
            int segmentEnd = segmentStart + segment.baseLength();
            if(segment instanceof RefSegment refSegment)
            {
                // Intersect the targeted range with this segment's span in probe-sequence space, then convert to offsets within the segment.
                int intersectionStart = max(segmentStart, targetedStart);
                int intersectionEnd = min(segmentEnd, targetedEnd);
                if(intersectionStart < intersectionEnd)
                {
                    targetedRegions.add(getSubregion(refSegment.region(), refSegment.orientation(),
                            intersectionStart - segmentStart, intersectionEnd - segmentStart));
                }
            }
            segmentStart = segmentEnd;
        }

        return targetedRegions;
    }
}
