package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_COVERAGE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionEndingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionStartingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.roundTowardsRegionEnds;

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

    // Calculates the maximum probe starting position such that the specified region overlaps the probe.
    public static int maxProbeStartOverlapping(final BaseRegion region)
    {
        return region.end();
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

    // Calculates the minimum probe starting position such that the probe overlaps an acceptable amount with the target.
    public static int minProbeStartPartiallyCovering(final BaseRegion target)
    {
        int minEnd = target.start() + PROBE_COVERAGE_MIN - 1;
        return probeRegionEndingAt(minEnd).start();
    }

    // Calculates the maximum probe ending position such that the probe overlaps an acceptable amount with the target.
    public static int maxProbeEndPartiallyCovering(final BaseRegion target)
    {
        int maxStart = target.end() - PROBE_COVERAGE_MIN + 1;
        return probeRegionStartingAt(maxStart).end();
    }

    // Calculates the best probe tiling of a region, respecting the probe overlap and coverage constraints.
    // Generally the probes are tiled such that they are centered on the region and equally spaced.
    // `probeBounds` indicates hard bounds that probes must be fully contained within. Otherwise, some extension is allowed.
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

        int maxUncovered = PROBE_LENGTH - PROBE_COVERAGE_MIN;
        // Hard bounds on the region in which probes can be placed.
        BaseRegion tilingBounds = new BaseRegion(
                max(probeBounds.start(), region.start() - maxUncovered),
                min(probeBounds.end(), region.end() + maxUncovered));

        double centre = (region.start() + region.end()) / 2.0;

        // Lower bound is number of probes which fit completely within the target region without overlap.
        int minProbes = region.baseLength() / PROBE_LENGTH;

        // Upper bound is maximally overlapped and maximally extending outside target region.
        // maxProbes * PROBE_LENGTH - (maxProbes-1) * PROBE_OVERLAP_MAX <= regionSize + startExtend + endExtend
        int maxProbes = (tilingBounds.baseLength() - PROBE_OVERLAP_MAX) / (PROBE_LENGTH - PROBE_OVERLAP_MAX);

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
            tilingStart = centre - tilingLength / 2;

            // Adjust the tiling alignment to adhere to the hard bounds.
            tilingStart = max(tilingStart, tilingBounds.start());
            double end = tilingStart + tilingLength / 2 - 1;
            if(end > tilingBounds.end())
            {
                tilingStart -= end - tilingBounds.end();
            }
        }

        double tilingStartCopy = tilingStart;
        List<Integer> probes = IntStream.range(0, probeCount)
                // Round towards ends of region to prefer grouping uncovered regions in the middle rather than missing a base at the end.
                .map(i -> (int) roundTowardsRegionEnds(tilingStartCopy + i * probeStartSpacing, region))
                .boxed().toList();
        return probes;
    }
}
