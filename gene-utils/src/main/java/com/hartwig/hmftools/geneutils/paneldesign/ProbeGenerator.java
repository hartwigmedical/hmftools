package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.canCoverRegionWithOneProbe;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeEndContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeEndCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.nextProbeStartPosition;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeStartingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.getBestScoringElement;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
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
    public ProbeGenerationResult coverRegion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        if(canCoverRegionWithOneProbe(target.region().baseRegion()))
        {
            return coverSmallRegion(target, criteria);
        }
        else
        {
            return coverLargeRegion(target, criteria);
        }
    }

    // Generates the best acceptable probes to cover an entire region larger than the probe size.
    private ProbeGenerationResult coverLargeRegion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        // Methodology:
        //   - For each position in the target, try to find the 1 best acceptable probe that covers it.
        //   - If an acceptable probe is found, advance the position to the next position after the probe.
        //   - If a position can't be covered by a probe, move to the next position.

        ChrBaseRegion targetBaseRegion = target.region();

        ProbeFactory probeFactory = new ProbeFactory(target);
        List<EvaluatedProbe> probes = new ArrayList<>();
        for(int position = targetBaseRegion.start(); position <= targetBaseRegion.end(); )
        {
            // Ensure the probe covers this position and doesn't overlap the previous probe.
            int minProbeStart = max(
                    minProbeStartContaining(position),
                    probes.isEmpty() ? 1 : nextProbeStartPosition(probes.get(probes.size() - 1).candidate().probeRegion().end()));
            Stream<CandidateProbe> candidates =
                    leftMovingLeftAlignedProbes(new BasePosition(targetBaseRegion.chromosome(), position), minProbeStart, probeFactory);
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
                targetBaseRegion.baseRegion(), probes.stream().map(probe -> probe.candidate().probeRegion().baseRegion()))
                .stream()
                .map(region ->
                        new RejectedRegion(ChrBaseRegion.from(targetBaseRegion.chromosome(), region), target, rejectionReason))
                .toList();

        return new ProbeGenerationResult(List.of(target), probes, rejectedRegions);
    }

    // Generates the 1 best acceptable probe that covers the specified smaller region.
    private ProbeGenerationResult coverSmallRegion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        return selectBestProbe(coverSmallRegionCandidates(target), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(target), List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(target),
                            Collections.emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(target, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> coverSmallRegionCandidates(final TargetRegion target)
    {
        ChrBaseRegion targetBaseRegion = target.region();

        if(!canCoverRegionWithOneProbe(targetBaseRegion.baseRegion()))
        {
            // This method is designed to find 1 probe covering the whole region, which is impossible in this case.
            throw new IllegalArgumentException("target must be no larger than probe");
        }

        int regionCentre = (targetBaseRegion.start() + targetBaseRegion.end()) / 2;
        BasePosition initialPosition = new BasePosition(targetBaseRegion.chromosome(), regionCentre);
        // Stop once the probes are too far from the target position to cover it.
        int minProbeStart = minProbeStartCovering(targetBaseRegion.baseRegion());
        int maxProbeEnd = maxProbeEndCovering(targetBaseRegion.baseRegion());
        ProbeFactory probeFactory = new ProbeFactory(target);
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, probeFactory);
    }

    // Generates the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult coverOneSubregion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        return selectBestProbe(coverOneSubregionCandidates(target), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(target), List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(target),
                            Collections.emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(target, rejectionReason)));
                });
    }

    private static Stream<CandidateProbe> coverOneSubregionCandidates(final TargetRegion target)
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
        ProbeFactory probeFactory = new ProbeFactory(target);
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd, probeFactory);
    }

    // Generate candidate probes which cover a position, starting from the position and moving outwards.
    public static Stream<CandidateProbe> coverPositionCandidates(final BasePosition position, final TargetMetadata metadata)
    {
        int minProbeStart = minProbeStartContaining(position.Position);
        int maxProbeEnd = maxProbeEndContaining(position.Position);
        ChrBaseRegion targetBaseRegion = ChrBaseRegion.from(position);
        TargetRegion target = new TargetRegion(targetBaseRegion, metadata);
        ProbeFactory probeFactory = new ProbeFactory(target);
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
