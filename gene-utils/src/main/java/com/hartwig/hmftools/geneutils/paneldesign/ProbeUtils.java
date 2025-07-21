package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_COVERAGE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionEndingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionStartingAt;

import java.util.List;
import java.util.stream.IntStream;

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

    // Calculates the best probe tiling of a region, respecting the probe overlap and coverage constraints.
    // `startExtend` and `endExtend` indicate if the probes may extend past the start and end of the region.
    public static List<BaseRegion> calculateIdealProbeTiling(final BaseRegion region, boolean startExtend, boolean endExtend)
    {
        // Lower bound is number of probes which fit completely within the target region without overlap.
        int minProbes = region.baseLength() / PROBE_LENGTH;

        // Upper bound is maximally overlapped and maximally extending outside target region.
        int maxUncovered = PROBE_LENGTH - PROBE_COVERAGE_MIN;
        int startUncoveredMax = startExtend ? maxUncovered : 0;
        int endUncoveredMax = endExtend ? maxUncovered : 0;
        // maxProbes * PROBE_LENGTH - (maxProbes-1) * PROBE_OVERLAP_MAX <= regionSize + startUncoveredMax + endUncoveredMax
        int maxProbes =
                (region.baseLength() + startUncoveredMax + endUncoveredMax - PROBE_OVERLAP_MAX) / (PROBE_LENGTH - PROBE_OVERLAP_MAX);

        // Select the probe count which minimises the number of overlapping/uncovered (extra) bases.
        // This is always minProbes or minProbes + 1, since minProbes minimally undershoots and minProbes+1 will cause positive extra bases.
        int probeCount = minProbes;
        int extra = minProbes * PROBE_LENGTH - region.baseLength();
        if(maxProbes > minProbes)
        {
            int nextExtra = (minProbes + 1) * PROBE_LENGTH - region.baseLength();
            // <= to prefer more probes which means more target coverage.
            if(abs(nextExtra) <= abs(extra))
            {
                probeCount = minProbes + 1;
                extra = nextExtra;
            }
        }

        if(probeCount <= 0)
        {
            // Can't even place a single probe with the given constraints.
            return emptyList();
        }
        else
        {
            // Calculate the start of the tiling and the space between each probe start position.
            // Then just tile the probes regularly from those parameters (since the constraints were checked earlier).
            double tilingStart;
            double probeStartSpacing;
            if(extra == 0)
            {
                // Perfect tiling.
                tilingStart = region.start();
                probeStartSpacing = PROBE_LENGTH;
            }
            else
            {
                // General case.
                probeStartSpacing = probeCount <= 1
                        ? 0
                        : max((region.baseLength() - PROBE_LENGTH) / (probeCount - 1), PROBE_LENGTH - PROBE_OVERLAP_MAX);
                double tilingLength = (probeCount - 1) * probeStartSpacing + PROBE_LENGTH;
                double regionCentreTrue = (region.start() + region.end()) / 2.0;
                // Pick tiling start position based on start and end bounds (if present).
                tilingStart =
                        startExtend ? (endExtend ? regionCentreTrue - tilingLength / 2 : region.end() - tilingLength + 1) : region.start();
            }

            return IntStream.range(0, probeCount)
                    .mapToObj(i -> probeRegionStartingAt((int) round(tilingStart + i * probeStartSpacing)))
                    .toList();
        }
    }
}
