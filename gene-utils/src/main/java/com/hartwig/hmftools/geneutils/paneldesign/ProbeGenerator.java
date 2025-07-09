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

    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final ProbeEvaluator probeEvaluator)
    {
        mProbeEvaluator = probeEvaluator;
    }

    public ProbeGenerationResult coverWholeRegion(final ChrBaseRegion region, final ProbeSourceInfo source)
    {
        // TODO
        return new ProbeGenerationResult();
    }

    // Find the 1 best probe that covers any acceptable part of the specified region.
    public ProbeGenerationResult bestProbeInRegion(final ChrBaseRegion region, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        return mProbeEvaluator.selectBestProbe(bestProbeInRegionCandidates(region, source), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                // TODO: rejection reason
                .orElseGet(() -> new ProbeGenerationResult(Collections.emptyList(), List.of(new RejectedRegion(region, source, null))));
    }

    private static Stream<CandidateProbe> bestProbeInRegionCandidates(final ChrBaseRegion region, final ProbeSourceInfo source)
    {
        // TODO
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public static Stream<CandidateProbe> coverPositionCandidates(final BasePosition position, final ProbeSourceInfo source)
    {
        ChrBaseRegion targetRegion = new ChrBaseRegion(position.Chromosome, position.Position, position.Position);
        // Generates probes whose centres are offset by these amounts: 0, +1, -1, +2, -2, +3, -3, ...
        return IntStream.iterate(0, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset))
                .mapToObj(offset ->
                {
                    int probeStart = position.Position + offset - PROBE_LENGTH / 2;
                    int probeEnd = probeStart + PROBE_LENGTH;
                    return new CandidateProbe(source, new ChrBaseRegion(position.Chromosome, probeStart, probeEnd), targetRegion);
                })
                // Negative offsets could go before
                .filter(probe -> probe.probeRegion().start() >= 1)
                .takeWhile(probe -> probeCoversTarget(probe.probeRegion().baseRegion(), targetRegion.baseRegion()))
                .peek(ProbeGenerator::logCandidateProbe);
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

    private static void logCandidateProbe(final CandidateProbe probe)
    {
        LOGGER.trace("Candidate probe: {}", probe);
    }
}
