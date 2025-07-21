package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeEndContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.outwardMovingOffsets;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentre;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentreStartOffset;

import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// TODO: unit test

// Utilities for creating candidate probes covering target regions.
public class CandidateProbeGenerator
{
    private final Map<String, Integer> mChromosomeLengths;

    public CandidateProbeGenerator(final Map<String, Integer> chromosomeLengths)
    {
        mChromosomeLengths = chromosomeLengths;
    }

    // Generates candidate probes covering any subregion within the target region.
    // Probes do not extend outside the target region.
    public Stream<CandidateProbe> coverOneSubregion(final ChrBaseRegion region, final CandidateProbeContext context)
    {
        if(region.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least 1 probe.
            throw new IllegalArgumentException("target must be larger than a probe");
        }

        BasePosition initialPosition = new BasePosition(region.chromosome(), regionCentre(region.baseRegion()));
        // Stop once the probes go outside the target region.
        int minProbeStart = region.start();
        int maxProbeEnd = region.end();
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, context);
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public Stream<CandidateProbe> coverPosition(final BasePosition position, final CandidateProbeContext context)
    {
        int minProbeStart = minProbeStartContaining(position.Position);
        int maxProbeEnd = maxProbeEndContaining(position.Position);
        return outwardMovingCenterAlignedProbes(position, minProbeStart, maxProbeEnd, context);
    }

    // Returns candidate probes shifting left and right with offsets: 0, 1, -1, 2, -2, 3, -3, ...
    // Probes are aligned to their centre position.
    // Probe bounds:
    //   - Can't start before start of chromosome
    //   - Can't start before `minProbeStart`
    //   - Can't end after `maxProbeEnd`
    //   - Can't end after end of chromosome
    // Useful because we prefer to select probes which are closest to the target position or centre of a region.
    public Stream<CandidateProbe> outwardMovingCenterAlignedProbes(final BasePosition initialPosition, int minProbeStart,
            int maxProbeEnd, final CandidateProbeContext context)
    {
        if(maxProbeEnd - minProbeStart + 1 < PROBE_LENGTH)
        {
            // Probably indicates a bug in the calling code.
            throw new IllegalArgumentException("minProbeStart and maxProbeEnd forbid all possible probes");
        }

        minProbeStart = max(minProbeStart, 1);
        maxProbeEnd = min(maxProbeEnd, mChromosomeLengths.get(initialPosition.Chromosome));

        // Must be consistent with probeRegionCenteredAt().
        int centreStartOffset = regionCentreStartOffset(PROBE_LENGTH);

        // minProbeStart = initialPosition + offset + centreStartOffset
        int minOffset = minProbeStart - initialPosition.Position - centreStartOffset;

        // maxProbeEnd = initialPosition + offset + centreStartOffset + PROBE_LENGTH - 1
        int maxOffset = maxProbeEnd - initialPosition.Position - centreStartOffset - PROBE_LENGTH + 1;

        return outwardMovingOffsets(minOffset, maxOffset)
                .mapToObj(offset -> context.createProbe(probeRegionCenteredAt(initialPosition.Position + offset)));
    }
}
