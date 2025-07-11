package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_OVERLAP_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.UNCOVERED_BASES_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.getBestScoringElement;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

public class ProbeGenerator
{
    public final ProbeEvaluator mProbeEvaluator;

    public ProbeGenerator(final ProbeEvaluator probeEvaluator)
    {
        mProbeEvaluator = probeEvaluator;
    }

    // TODO: check chromosome lengths and reject probes outside chromosome

    // Generates the best acceptable probes to cover an entire region. The probes may overlap and extend outside the target region.
    public ProbeGenerationResult coverRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        if(canCoverRegionWith1Probe(target))
        {
            return coverSmallRegion(target, source, criteria);
        }
        else
        {
            return coverLargeRegion(target, source, criteria);
        }
    }

    // Generates the best acceptable probes to cover an entire region larger than the probe size.
    private ProbeGenerationResult coverLargeRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        // Methodology:
        //   - For each position in the target, try to find the 1 best acceptable probe that covers it.
        //   - If an acceptable probe is found, advance the position to the next position after the probe.
        //   - If a position can't be covered by a probe, move to the next position.

        ProbeFactory probeFactory = new ProbeFactory(source, target);
        List<EvaluatedProbe> probes = new ArrayList<>();
        for(int position = target.start(); position <= target.end(); )
        {
            // Ensure the probe covers this position and doesn't overlap the previous probe.
            int minProbeStart = max(
                    minProbeStartContaining(position),
                    probes.isEmpty() ? 1 : nextProbeStartPosition(probes.get(probes.size() - 1).candidate().probeRegion().end()));
            Stream<CandidateProbe> candidates =
                    leftMovingLeftAlignedProbes(new BasePosition(target.chromosome(), position), minProbeStart, probeFactory);
            Optional<EvaluatedProbe> bestCandidate = selectBestProbe(candidates, criteria);

            if(bestCandidate.isPresent())
            {
                probes.add(bestCandidate.get());
                position = bestCandidate.get().candidate().probeRegion().end() + 1;
            }
            else
            {
                ++position;
            }
        }

        // Compute rejected regions based on what has been covered by the probes.
        String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
        List<RejectedRegion> rejectedRegions = computeUncoveredRegions(
                target.baseRegion(), probes.stream().map(probe -> probe.candidate().probeRegion().baseRegion()))
                .stream()
                .map(region ->
                        new RejectedRegion(new ChrBaseRegion(target.chromosome(), region.start(), region.end()), source, rejectionReason))
                .toList();

        return new ProbeGenerationResult(probes, rejectedRegions);
    }

    // Generates the 1 best acceptable probe that covers the specified smaller region.
    private ProbeGenerationResult coverSmallRegion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        return selectBestProbe(coverSmallRegionCandidates(target, source), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            Collections.emptyList(),
                            List.of(new RejectedRegion(target, source, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> coverSmallRegionCandidates(final ChrBaseRegion target, final ProbeSourceInfo source)
    {
        if(!canCoverRegionWith1Probe(target))
        {
            // This method is designed to find 1 probe covering the whole region, which is impossible in this case.
            throw new IllegalArgumentException("target must be no larger than probe");
        }

        int regionCentre = (target.start() + target.end()) / 2;
        BasePosition initialPosition = new BasePosition(target.chromosome(), regionCentre);
        // Stop once the probes are too far from the target position to cover it.
        int minProbeStart = minProbeStartCovering(target.baseRegion());
        int maxProbeEnd = maxProbeEndCovering(target.baseRegion());
        ProbeFactory probeFactory = new ProbeFactory(source, target);
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, probeFactory);
    }

    // Generates the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult coverOneSubregion(final ChrBaseRegion target, final ProbeSourceInfo source,
            final ProbeSelectCriteria criteria)
    {
        return selectBestProbe(coverOneSubregionCandidates(target, source), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            Collections.emptyList(),
                            List.of(new RejectedRegion(target, source, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> coverOneSubregionCandidates(final ChrBaseRegion target, final ProbeSourceInfo source)
    {
        if(target.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least 1 probe.
            throw new IllegalArgumentException("target must be larger than a probe");
        }

        int regionCentre = (target.start() + target.end()) / 2;
        BasePosition initialPosition = new BasePosition(target.chromosome(), regionCentre);
        // Stop once the probes go outside the target region.
        int minProbeStart = target.start();
        int maxProbeEnd = target.end();
        ProbeFactory probeFactory = new ProbeFactory(source, target);
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, probeFactory);
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public static Stream<CandidateProbe> coverPositionCandidates(final BasePosition position, final ProbeSourceInfo source)
    {
        int minProbeStart = minProbeStartContaining(position.Position);
        int maxProbeEnd = maxProbeEndContaining(position.Position);
        ChrBaseRegion targetRegion = new ChrBaseRegion(position.Chromosome, position.Position, position.Position);
        ProbeFactory probeFactory = new ProbeFactory(source, targetRegion);
        return outwardMovingCenterAlignedProbes(position, minProbeStart, maxProbeEnd, probeFactory);
    }

    // Returns candidate probes shifting left and right with offsets: 0, 1, -1, 2, -2, 3, -3, ...
    // Probes are aligned to their centre position.
    // Probe bounds:
    //   - Can't start before start of chromosome
    //   - Can't start before `minProbeStart`
    //   - Can't end after `maxProbeEnd`
    // TODO: don't extend past end of chromosome
    // Useful because we prefer to select probes which are closest to the target position or centre of a region.
    private static Stream<CandidateProbe> outwardMovingCenterAlignedProbes(final BasePosition initialPosition, int minProbeStart,
            int maxProbeEnd, final ProbeFactory factory)
    {
        minProbeStart = max(1, minProbeStart);

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
    // TODO: don't extend past end of chromosome
    private static Stream<CandidateProbe> leftMovingLeftAlignedProbes(final BasePosition position, int minProbeStart,
            final ProbeFactory factory)
    {
        minProbeStart = max(1, minProbeStart);
        // minProbeStart <= position + offset
        int minOffset = minProbeStart - position.Position;
        int maxOffset = 0;
        return IntStream.iterate(maxOffset, offset -> offset >= minOffset, offset -> offset - 1)
                .mapToObj(offset -> probeStartingAt(position.Chromosome, position.Position + offset, factory));
    }

    private static boolean canCoverRegionWith1Probe(final ChrBaseRegion target)
    {
        return target.baseLength() - UNCOVERED_BASES_MAX <= PROBE_LENGTH;
    }

    private static int minProbeStartCovering(final BaseRegion target)
    {
        // Need to take care of the case where target is smaller than a probe: don't allow the probe to end before the target start.
        int minProbeEnd = max(target.start(), target.end() - UNCOVERED_BASES_MAX);
        return minProbeEnd - PROBE_LENGTH + 1;
    }

    private static int maxProbeEndCovering(final BaseRegion target)
    {
        // Need to take care of the case where target is smaller than a probe: don't allow the probe to start after the target end.
        int maxProbeStart = min(target.end(), target.start() + UNCOVERED_BASES_MAX);
        return maxProbeStart + PROBE_LENGTH - 1;
    }

    // Calculates the minimum probe starting position such that the specified position is contained within the probe.
    private static int minProbeStartContaining(int targetPosition)
    {
        // start + PROBE_LENGTH - 1 >= targetPosition
        return targetPosition - PROBE_LENGTH + 1;
    }

    // Calculates the maximum probe ending position such that the specified position is contained within the probe.
    private static int maxProbeEndContaining(int targetPosition)
    {
        // end - PROBE_LENGTH + 1 <= targetPosition
        return targetPosition + PROBE_LENGTH - 1;
    }

    // Calculates the next position a probe may start at, respecting the probe overlap constraint.
    private static int nextProbeStartPosition(int prevProbeEnd)
    {
        // prevProbeEnd - nextStart + 1 < PROBE_OVERLAP_MAX
        return prevProbeEnd - PROBE_OVERLAP_MAX + 1;
    }

    // Encapsulates the constant fields so we don't have to pass them around everywhere when creating candidates.
    private record ProbeFactory(
            ProbeSourceInfo source,
            ChrBaseRegion targetRegion
    )
    {
        public CandidateProbe create(final ChrBaseRegion probeRegion)
        {
            return new CandidateProbe(source, probeRegion, targetRegion);
        }
    }

    private static CandidateProbe probeCenteredAt(final String chromosome, int centrePosition, final ProbeFactory factory)
    {
        return probeStartingAt(chromosome, centrePosition - PROBE_LENGTH / 2, factory);
    }

    private static CandidateProbe probeStartingAt(final String chromosome, int startPosition, final ProbeFactory factory)
    {
        return factory.create(new ChrBaseRegion(chromosome, startPosition, startPosition + PROBE_LENGTH - 1));
    }

    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there are no acceptable probes.
    public Optional<EvaluatedProbe> selectBestProbe(Stream<CandidateProbe> probes, final ProbeSelectCriteria criteria)
    {
        Stream<EvaluatedProbe> acceptableProbes = probes
                .map(probe -> mProbeEvaluator.evaluateCandidate(probe, criteria.eval()))
                .filter(EvaluatedProbe::accepted);
        return switch(criteria.strategy())
        {
            case FIRST_ACCEPTABLE -> acceptableProbes.findFirst();
            // Stop if a probe with quality=1 is found.
            case MAX_QUALITY -> getBestScoringElement(acceptableProbes,
                    EvaluatedProbe::qualityScore,
                    quality -> Doubles.greaterOrEqual(quality, 1.0), true);
            // Stop if a probe with gc=gcTarget is found.
            case BEST_GC -> getBestScoringElement(acceptableProbes,
                    probe -> abs(probe.gcContent() - criteria.eval().gcContentTarget()),
                    gc -> Doubles.equal(gc, criteria.eval().gcContentTarget()), false);
        };
    }
}
