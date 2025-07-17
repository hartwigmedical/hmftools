package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeEndContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeStartingAt;

import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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
    public Stream<CandidateProbe> coverOneSubregion(final TargetRegion target)
    {
        ChrBaseRegion targetBaseRegion = target.region();

        if(targetBaseRegion.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least 1 probe.
            throw new IllegalArgumentException("target must be larger than a probe");
        }

        int regionCentre = (targetBaseRegion.start() + targetBaseRegion.end()) / 2;
        BasePosition initialPosition = new BasePosition(targetBaseRegion.chromosome(), regionCentre);
        // Stop once the probes go outside the target region.
        int minProbeStart = targetBaseRegion.start();
        int maxProbeEnd = targetBaseRegion.end();
        CandidateProbeFactory candidateFactory = new CandidateProbeFactory(target);
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, candidateFactory);
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public Stream<CandidateProbe> coverPosition(final BasePosition position, final TargetMetadata metadata)
    {
        int minProbeStart = minProbeStartContaining(position.Position);
        int maxProbeEnd = maxProbeEndContaining(position.Position);
        ChrBaseRegion targetBaseRegion = ChrBaseRegion.from(position);
        TargetRegion target = new TargetRegion(targetBaseRegion, metadata);
        CandidateProbeFactory candidateFactory = new CandidateProbeFactory(target);
        return outwardMovingCenterAlignedProbes(position, minProbeStart, maxProbeEnd, candidateFactory);
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
            int maxProbeEnd, final CandidateProbeFactory factory)
    {
        if(maxProbeEnd - minProbeStart < PROBE_LENGTH)
        {
            // Probably indicates a bug in the calling code.
            throw new IllegalArgumentException("minProbeStart and maxProbeEnd forbid all possible probes");
        }

        minProbeStart = max(minProbeStart, 1);
        maxProbeEnd = min(maxProbeEnd, mChromosomeLengths.get(initialPosition.Chromosome));

        // minProbeStart = initialPosition + offset - PROBE_LENGTH/2
        // minProbeStart - initialPosition + PROBE_LENGTH/2 = offset
        int minOffset = minProbeStart - initialPosition.Position + PROBE_LENGTH / 2;

        // maxProbeEnd = initialPosition + offset - PROBE_LENGTH/2 + PROBE_LENGTH - 1;
        // maxProbeEnd - initialPosition + PROBE_LENGTH/2 - PROBE_LENGTH + 1 = offset
        int maxOffset = maxProbeEnd - initialPosition.Position + PROBE_LENGTH / 2 - PROBE_LENGTH + 1;

        return IntStream.iterate(0, absOffset -> -absOffset >= minOffset || absOffset <= maxOffset, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset))
                .filter(offset -> offset >= minOffset && offset <= maxOffset)
                .mapToObj(offset -> probeCenteredAt(initialPosition.Chromosome, initialPosition.Position + offset, factory));
    }

    // Generates probes shifting with offsets: 0, -1, -2, -3, ...
    // Probes are aligned to the start position.
    // Probe bounds:
    //   - Can't start before start of chromosome
    //   - Can't start before `minProbeStart`
    //   - Can't end after end of chromosome
    public Stream<CandidateProbe> leftMovingLeftAlignedProbes(final BasePosition initialPosition, int minProbeStart,
            final CandidateProbeFactory factory)
    {
        if(initialPosition.Position < minProbeStart)
        {
            // Probably indicates a bug in the calling code.
            throw new IllegalArgumentException("minProbeStart forbids all possible probes");
        }

        minProbeStart = max(minProbeStart, 1);
        int maxProbeEnd = min(initialPosition.Position + PROBE_LENGTH - 1, mChromosomeLengths.get(initialPosition.Chromosome));

        // minProbeStart = initialPosition + offset
        int minOffset = minProbeStart - initialPosition.Position;

        // maxProbeEnd = initialPosition + offset + PROBE_LENGTH - 1
        int maxOffset = maxProbeEnd - initialPosition.Position - PROBE_LENGTH + 1;

        return IntStream.iterate(maxOffset, offset -> offset >= minOffset, offset -> offset - 1)
                .mapToObj(offset -> probeStartingAt(initialPosition.Chromosome, initialPosition.Position + offset, factory));
    }
}
