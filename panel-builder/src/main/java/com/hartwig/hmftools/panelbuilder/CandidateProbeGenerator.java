package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentre;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentreStartOffset;
import static com.hartwig.hmftools.panelbuilder.Utils.outwardMovingOffsets;

import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Functionality for creating candidate probes covering target regions.
public class CandidateProbeGenerator
{
    private final Map<String, Integer> mChromosomeLengths;

    public CandidateProbeGenerator(final Map<String, Integer> chromosomeLengths)
    {
        mChromosomeLengths = chromosomeLengths;
    }

    // Generates candidate probes covering any subregion within the target region.
    // Probes do not extend outside the target region.
    public Stream<Probe> coverOneSubregion(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        if(region.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least one probe.
            throw new IllegalArgumentException("region must be larger than a probe");
        }

        BasePosition initialPosition = regionCentre(region);
        // Stop once the probes go outside the target region.
        int minProbeStart = region.start();
        int maxProbeEnd = region.end();
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, metadata);
    }

    // Returns candidate probes shifting left and right with offsets: 0, 1, -1, 2, -2, 3, -3, ...
    // Probes are aligned to their centre position.
    // Probe bounds:
    //   - Can't start before start of chromosome
    //   - Can't start before `minProbeStart`
    //   - Can't end after `maxProbeEnd`
    //   - Can't end after end of chromosome
    // Useful because we prefer to select probes which are closest to the target position or centre of a region.
    private Stream<Probe> outwardMovingCenterAlignedProbes(final BasePosition initialPosition, int minProbeStart, int maxProbeEnd,
            final TargetMetadata metadata)
    {
        if(maxProbeEnd - minProbeStart + 1 < PROBE_LENGTH)
        {
            // Probably indicates a bug in the calling code.
            throw new IllegalArgumentException("minProbeStart and maxProbeEnd forbid all possible probes");
        }

        // Must be consistent with probeRegionCenteredAt().
        int centreStartOffset = regionCentreStartOffset(PROBE_LENGTH);

        // minProbeStart = initialPosition + offset + centreStartOffset
        int minOffset = minProbeStart - initialPosition.Position - centreStartOffset;

        // maxProbeEnd = initialPosition + offset + centreStartOffset + PROBE_LENGTH - 1
        int maxOffset = maxProbeEnd - initialPosition.Position - centreStartOffset - PROBE_LENGTH + 1;

        return outwardMovingOffsets(minOffset, maxOffset)
                .mapToObj(offset ->
                {
                    SequenceDefinition definition = SequenceDefinition.singleRegion(
                            probeRegionCenteredAt(initialPosition.Chromosome, initialPosition.Position + offset));
                    return new Probe(definition, metadata);
                });
    }

    // Generates all probes overlapping a region, in order from left to right.
    // Use with care. Generally you would want to choose probes with a more careful approach.
    public Stream<Probe> allOverlapping(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        int minProbeStart = max(minProbeStartOverlapping(region.baseRegion()), 1);
        int maxProbeEnd = min(maxProbeEndOverlapping(region.baseRegion()), mChromosomeLengths.get(region.chromosome()));
        int maxProbeStart = probeRegionEndingAt(maxProbeEnd).start();
        return IntStream.rangeClosed(minProbeStart, maxProbeStart)
                .mapToObj(start ->
                {
                    SequenceDefinition definition = SequenceDefinition.singleRegion(probeRegionStartingAt(region.chromosome(), start));
                    return new Probe(definition, metadata);
                });
    }
}
