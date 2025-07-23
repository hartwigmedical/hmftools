package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentreFloat;
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

    // Calculates the best probe tiling of a region.
    // Objectives:
    //   - Cover the whole region;
    //   - Centre the tiling on the region;
    //   - Equally spaced probes;
    //   - Probes cannot extend outside `probeBounds`.
    // Returns the start positions of the probes.
    public static List<Integer> calculateOptimalProbeTiling(final BaseRegion region, final BaseRegion probeBounds)
    {
        if(!region.hasValidPositions() || !probeBounds.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(!probeBounds.containsRegion(region))
        {
            // Probably a bug in the caller.
            throw new IllegalArgumentException("probeBounds forbids all possible tilings");
        }

        int maxUncovered = PROBE_LENGTH - 1;
        // Hard bounds on the region in which probes can be placed.
        BaseRegion tilingBounds = new BaseRegion(
                max(probeBounds.start(), region.start() - maxUncovered),
                min(probeBounds.end(), region.end() + maxUncovered));

        double centre = regionCentreFloat(region);

        // Lower bound is number of probes which completely cover the region.
        int minProbes = (int) ceil((double) region.baseLength() / PROBE_LENGTH);

        // Upper bound is maximally overlapped and maximally extending outside the target region.
        int maxOverlap = PROBE_LENGTH - 1;
        // maxProbes * PROBE_LENGTH - (maxProbes-1) * maxOverlap <= regionSize + startExtend + endExtend
        int maxProbes = (tilingBounds.baseLength() - maxOverlap) / (PROBE_LENGTH - maxOverlap);

        if(minProbes > maxProbes)
        {
            // No tiling is possible given the constraints.
            return emptyList();
        }

        // The optimal number of probes is always the minimum possible, since we guaranteed the minimum covers the whole region, and adding
        // more probes will only increase overlap or extension for no gain.
        int probeCount = minProbes;
        // How many probe bases "left over" from covering the region?
        int extra = minProbes * PROBE_LENGTH - region.baseLength();

        if(probeCount <= 0)
        {
            // Can't even place a single probe with the given constraints.
            return emptyList();
        }

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
            // TODO: parameter to balance overlap and extension
            probeStartSpacing = probeCount <= 1
                    ? 0.0
                    : max((region.baseLength() - PROBE_LENGTH) / (probeCount - 1.0), PROBE_LENGTH - maxOverlap);
            double tilingLength = (probeCount - 1) * probeStartSpacing + PROBE_LENGTH;
            tilingStart = centre - tilingLength / 2;

            // Adjust the tiling alignment to adhere to the hard bounds.
            tilingStart = max(tilingStart, tilingBounds.start());
            double end = tilingStart + tilingLength - 1;
            if(end > tilingBounds.end())
            {
                tilingStart -= end - tilingBounds.end();
            }
        }

        double tilingStartCopy = tilingStart;
        List<Integer> probes = IntStream.range(0, probeCount)
                .map(i -> (int) round(tilingStartCopy + i * probeStartSpacing))
                .boxed().toList();
        return probes;
    }
}
