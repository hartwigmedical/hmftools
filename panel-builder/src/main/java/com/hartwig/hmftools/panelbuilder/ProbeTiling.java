package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_OVERLAP_EXTENSION_BALANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REGION_UNCOVERED_MAX;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentreFloat;

import java.util.List;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.region.BaseRegion;

// Pure positional maths for evenly tiling probes over a 1-D position range. Coordinate-space agnostic, so shared by DNA and RNA tiling.
public class ProbeTiling
{
    // Calculates the best probe tiling of a region.
    // Objectives:
    //   - Cover most of the region, with possibly a few bases on the edge uncovered;
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

        int maxProbeOffTarget = PROBE_LENGTH - 1;
        // Hard bounds on the region in which probes can be placed.
        BaseRegion tilingBounds = new BaseRegion(
                max(probeBounds.start(), region.start() - maxProbeOffTarget),
                min(probeBounds.end(), region.end() + maxProbeOffTarget));

        double centre = regionCentreFloat(region);

        // Lower bound is number of probes which completely cover the region, possibly excluding the max allowed uncovered bases.
        int minProbes = (int) ceil(max(1.0, region.baseLength() - REGION_UNCOVERED_MAX) / PROBE_LENGTH);

        // Upper bound is maximally overlapped and maximally extending outside the target region.
        int maxProbeOverlap = PROBE_LENGTH - 1;
        // maxProbes * PROBE_LENGTH - (maxProbes-1) * maxProbeOverlap <= regionSize + startExtend + endExtend
        int maxProbes = (tilingBounds.baseLength() - maxProbeOverlap) / (PROBE_LENGTH - maxProbeOverlap);

        if(minProbes > maxProbes)
        {
            // No tiling is possible given the constraints.
            return emptyList();
        }

        // The optimal number of probes is always the minimum possible, since we guaranteed the minimum acceptably covers the region,
        // and adding more probes will only increase overlap or extension for no gain.
        int probeCount = minProbes;
        // How many probe bases "left over" from covering the region?
        int extra = probeCount * PROBE_LENGTH - region.baseLength();

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

            // If extra is negative (i.e. some bases uncovered) then reduce the desired tiling region.
            // If extra is positive (i.e. some overlap) then distribute that between overlap and extension based on the constant.
            double tilingSpace = min(
                    region.baseLength() + min(0, extra) + max(0, extra * PROBE_OVERLAP_EXTENSION_BALANCE),
                    tilingBounds.baseLength());
            probeStartSpacing = probeCount <= 1
                    ? 0.0
                    : max(
                            (tilingSpace - PROBE_LENGTH) / (probeCount - 1.0),
                            PROBE_LENGTH - maxProbeOverlap);
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

    // Unified tiling entry point. With neither edge pinned this is the ordinary centred tiling (the DNA path); pinning an edge flush to an
    // exon boundary is the RNA generalisation for splice-junction coverage.
    public static List<Integer> calculateProbeTiling(final BaseRegion region, final BaseRegion probeBounds, boolean pinStart, boolean pinEnd)
    {
        if(!pinStart && !pinEnd)
        {
            return calculateOptimalProbeTiling(region, probeBounds);
        }
        // A pinned edge is flush to the boundary, so probes stay contained within the region (no extension past the pinned edge).
        return calculateContainedTiling(region, pinStart, pinEnd);
    }

    // Tiles probes strictly contained within the region, optionally pinning the outermost probes flush to the region edges (unlike the centred
    // tiling, probes never extend outside). Returns the probe start positions, empty if the region is shorter than a probe.
    public static List<Integer> calculateContainedTiling(final BaseRegion region, boolean pinStart, boolean pinEnd)
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }

        int length = region.baseLength();
        if(length < PROBE_LENGTH)
        {
            return emptyList();
        }

        // The uncovered budget is for the whole region, not per edge, and a pinned edge must be flush - so it applies only when an edge is
        // unpinned.
        int allowedUncovered = (pinStart && pinEnd) ? 0 : REGION_UNCOVERED_MAX;
        int probeCount = (int) ceil(max(1, length - allowedUncovered) / (double) PROBE_LENGTH);
        // How much of the region the probes span; the remainder (if any) is the allowed uncovered edge slack.
        int coveredLength = min(probeCount * PROBE_LENGTH, length);

        // Offset of the first probe within the region: flush left if the start is pinned, flush right if only the end is pinned, else
        // centred (splitting the slack between the two unpinned edges).
        int blockStart;
        if(pinStart)
        {
            blockStart = 0;
        }
        else if(pinEnd)
        {
            blockStart = length - coveredLength;
        }
        else
        {
            blockStart = (int) round((length - coveredLength) / 2.0);
        }

        double spacing = probeCount <= 1 ? 0.0 : (coveredLength - PROBE_LENGTH) / (double) (probeCount - 1);
        return IntStream.range(0, probeCount)
                .map(i -> region.start() + blockStart + (int) round(i * spacing))
                .boxed().toList();
    }
}
