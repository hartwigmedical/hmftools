package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.UNCOVERED_BASES_MAX;

import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ProbeGenerator
{
    public final ProbeEvaluator mProbeEvaluator;

    // TODO: needed?
    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final ProbeEvaluator probeEvaluator)
    {
        mProbeEvaluator = probeEvaluator;
    }

    public ProbeGenerationResult coverWholeRegion(final ChrBaseRegion region, final ProbeSourceInfo source)
    {
        // TODO
        return new ProbeGenerationResult(Collections.emptyList(), List.of(new RejectedRegion(region, source, "unimplemented algorithm")));
    }

    // Find the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult bestProbeInRegion(final ChrBaseRegion region, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        return mProbeEvaluator.selectBestProbe(bestProbeInRegionCandidates(region, source), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            Collections.emptyList(),
                            List.of(new RejectedRegion(region, source, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> bestProbeInRegionCandidates(final ChrBaseRegion region, final ProbeSourceInfo source)
    {
        if(region.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least 1 probe.
            throw new IllegalArgumentException("region must be larger than a probe");
        }

        int regionCentre = (region.start() + region.end()) / 2;
        return outwardMovingOffsets()
                .mapToObj(offset -> probeCenteredAt(region.chromosome(), regionCentre + offset, region, source))
                // Negative offsets could go before the start of the chromosome.
                .filter(probe -> probe.probeRegion().start() >= 1)
                // Stop once the probes go outside the target region.
                .takeWhile(probe -> region.containsRegion(probe.probeRegion()));
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public static Stream<CandidateProbe> coverPositionCandidates(final BasePosition position, final ProbeSourceInfo source)
    {
        ChrBaseRegion targetRegion = new ChrBaseRegion(position.Chromosome, position.Position, position.Position);
        return outwardMovingOffsets()
                .mapToObj(offset -> probeCenteredAt(position.Chromosome, position.Position + offset, targetRegion, source))
                // Negative offsets could go before the start of the chromosome.
                .filter(probe -> probe.probeRegion().start() >= 1)
                // Stop once the probes are too far from the target position to cover it.
                .takeWhile(probe -> probeCoversTarget(probe.probeRegion().baseRegion(), targetRegion.baseRegion()));
    }

    // Returns the sequence: 0, 1, -1, 2, -2, 3, -3, ...
    private static IntStream outwardMovingOffsets()
    {
        return IntStream.iterate(0, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset));
    }

    // Checks if a probe acceptably covers a target region.
    private static boolean probeCoversTarget(final BaseRegion probe, final BaseRegion target)
    {
        int uncoveredLeft = max(0, probe.start() - target.start());
        int uncoveredRight = max(0, target.end() - probe.end());
        int uncovered = uncoveredLeft + uncoveredRight;
        // First clause covers the case where target is smaller than probe.
        return uncovered < target.baseLength() && uncovered <= UNCOVERED_BASES_MAX;
    }

    private static CandidateProbe probeCenteredAt(final String chromosome, int centrePosition, final ChrBaseRegion targetRegion,
            final ProbeSourceInfo source)
    {
        return probeStartingAt(chromosome, centrePosition - PROBE_LENGTH / 2, targetRegion, source);
    }

    private static CandidateProbe probeStartingAt(final String chromosome, int startPosition, final ChrBaseRegion targetRegion,
            final ProbeSourceInfo source)
    {
        return new CandidateProbe(source, new ChrBaseRegion(chromosome, startPosition, startPosition + PROBE_LENGTH), targetRegion);
    }
}
