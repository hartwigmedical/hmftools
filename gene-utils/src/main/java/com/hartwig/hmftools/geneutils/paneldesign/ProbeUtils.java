package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_COVERAGE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentre;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionEndingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionStartingAt;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

// TODO: unit test

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

    // minProbeStartContaining() and maxProbeEndContaining() in one call.
    public static BaseRegion probeBoundsContaining(int targetPosition)
    {
        return new BaseRegion(minProbeStartContaining(targetPosition), maxProbeEndContaining(targetPosition));
    }

    // Calculates the next position a probe may start at, respecting the probe overlap constraint.
    public static int nextProbeStartPosition(int prevProbeEnd)
    {
        return prevProbeEnd - PROBE_OVERLAP_MAX + 1;
    }

    // Calculates the previous position a probe may end at, respecting the probe overlap constraint.
    public static int prevProbeEndPosition(int nextProbeStart)
    {
        return nextProbeStart + PROBE_OVERLAP_MAX - 1;
    }

    public static List<BaseRegion> calculateIdealProbeTiling(final BaseRegion region)
    {
        // TODO: this is not great

        // There are two probe patterns:
        //   - Odd number of probes, where the middle probe is centered on the region centre
        //   - Even number of probes, where the 2 middle probes are adjacent to the region centre
        // For now, only consider the most tightly tiled odd pattern.

        List<BaseRegion> result = new ArrayList<>();

        // Add centre probe.
        result.add(probeRegionCenteredAt(regionCentre(region)));
        // Try to add probes on to ends.
        while(true)
        {
            boolean change = false;
            BaseRegion first = result.get(0);
            if(first.start() > region.start())
            {
                // Try to add a probe to the left.
                BaseRegion prev = probeRegionEndingAt(prevProbeEndPosition(first.start()));
                int coverage = new BaseRegion(max(prev.start(), region.start()), prev.end()).baseLength();
                if(coverage >= PROBE_COVERAGE_MIN)
                {
                    result.add(0, prev);
                    change = true;
                }
            }
            BaseRegion last = result.get(result.size() - 1);
            if(last.end() < region.end())
            {
                // Try to add a probe to the right.
                BaseRegion next = probeRegionStartingAt(nextProbeStartPosition(last.end()));
                int coverage = new BaseRegion(next.start(), min(next.end(), region.end())).baseLength();
                if(coverage >= PROBE_COVERAGE_MIN)
                {
                    result.add(next);
                    change = true;
                }
            }
            if(!change)
            {
                break;
            }
        }

        // TODO: expand probes to fill region

        return result;
    }
}
