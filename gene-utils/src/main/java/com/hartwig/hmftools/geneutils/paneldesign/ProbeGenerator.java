package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.UNCOVERED_BASES_MAX;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
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

    public ProbeGenerationResult coverWholeRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeEvalCriteria criteria)
    {
        // TODO: probe eval/select criteria?
        // TODO: need to have different algorithms here?
        if(canCoverRegionWith1Probe(target))
        {
            // TODO: select strategy?
            return bestProbeCoveringRegion(target, source, new ProbeSelectCriteria(criteria, ProbeSelectStrategy.MAX_QUALITY));
        }
        else
        {
            // TODO
            return new ProbeGenerationResult(Collections.emptyList(), List.of(new RejectedRegion(target, source, "unimplemented algorithm")));
        }
    }

    // Find the 1 best acceptable probe that covers the specified smaller region.
    private ProbeGenerationResult bestProbeCoveringRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        Optional<EvaluatedProbe> bestCandidate =
                mProbeEvaluator.selectBestProbe(bestProbeCoveringRegionCandidates(target, source), criteria);
        return bestCandidate
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            Collections.emptyList(),
                            List.of(new RejectedRegion(target, source, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> bestProbeCoveringRegionCandidates(final ChrBaseRegion target, final ProbeSourceInfo source)
    {
        if(!canCoverRegionWith1Probe(target))
        {
            // This method is designed to find 1 probe covering the whole region, which is impossible in this case.
            throw new IllegalArgumentException("target must be no larger than probe");
        }

        int regionCentre = (target.start() + target.end()) / 2;
        return outwardMovingCenteredProbes(new BasePosition(target.chromosome(), regionCentre), target, source)
                // Stop once the probes are too far from the target position to cover it.
                .takeWhile(probe -> probeCoversTarget(probe.probeRegion().baseRegion(), target.baseRegion()));
    }

    // Find the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult bestProbeInRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        return mProbeEvaluator.selectBestProbe(bestProbeInRegionCandidates(target, source), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            Collections.emptyList(),
                            List.of(new RejectedRegion(target, source, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> bestProbeInRegionCandidates(final ChrBaseRegion target, final ProbeSourceInfo source)
    {
        if(target.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least 1 probe.
            throw new IllegalArgumentException("target must be larger than a probe");
        }

        int regionCentre = (target.start() + target.end()) / 2;
        return outwardMovingCenteredProbes(new BasePosition(target.chromosome(), regionCentre), target, source)
                // Stop once the probes go outside the target region.
                .takeWhile(probe -> target.containsRegion(probe.probeRegion()));
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public static Stream<CandidateProbe> coverPositionCandidates(final BasePosition position, final ProbeSourceInfo source)
    {
        ChrBaseRegion targetRegion = new ChrBaseRegion(position.Chromosome, position.Position, position.Position);
        return outwardMovingCenteredProbes(position, targetRegion, source)
                // Stop once the probes are too far from the target position to cover it.
                .takeWhile(probe -> probeCoversTarget(probe.probeRegion().baseRegion(), targetRegion.baseRegion()));
    }

    // Returns candidate probes starting centered on a position and then shifting left and right with offsets: 0, 1, -1, 2, -2, 3, -3, ...
    // Useful because we prefer to select probes which are closest to the target position or centre of a region.
    private static Stream<CandidateProbe> outwardMovingCenteredProbes(final BasePosition position, final ChrBaseRegion targetRegion,
            final ProbeSourceInfo source)
    {
        return IntStream.iterate(0, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset))
                .mapToObj(offset -> probeCenteredAt(position.Chromosome, position.Position + offset, targetRegion, source))
                // Negative offsets could go before the start of the chromosome.
                .filter(probe -> probe.probeRegion().start() >= 1);
    }

    private static boolean canCoverRegionWith1Probe(final ChrBaseRegion target)
    {
        return target.baseLength() - UNCOVERED_BASES_MAX <= PROBE_LENGTH;
    }

    // Checks if a probe acceptably covers a target region.
    private static boolean probeCoversTarget(final BaseRegion probe, final BaseRegion target)
    {
        // TODO: only ensure >= 1 base of the probe overlaps the target. Want to be stricter?
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
