package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.UNCOVERED_BASES_MAX;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.BiPredicate;
import java.util.function.DoublePredicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

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
        //   - If a position can't be covered by a probe, mark it as rejected and move to the next position.

        List<EvaluatedProbe> probes = new ArrayList<>();
        List<Integer> rejectedPositions = new ArrayList<>();
        for(int windowStart = target.start(); windowStart <= target.end(); )
        {
            // TODO: probe may significantly overlap last iteration. Introduce max overlap constraint?
            // TODO: reuse other code?
            int windowStartCopy = windowStart;
            Stream<CandidateProbe> candidates = IntStream.iterate(0, offset -> offset - 1)
                    .mapToObj(offset -> probeStartingAt(target.chromosome(), windowStartCopy + offset, target, source))
                    // Negative offsets could go before the start of the chromosome.
                    .filter(probe -> probe.probeRegion().start() >= 1)
                    .takeWhile(probe -> probe.probeRegion().containsPosition(windowStartCopy));
            Optional<EvaluatedProbe> bestCandidate = selectBestProbe(candidates, criteria);

            if(bestCandidate.isPresent())
            {
                probes.add(bestCandidate.get());
                windowStart = bestCandidate.get().candidate().probeRegion().end() + 1;
            }
            else
            {
                // We have checked up to windowStart-1. And in subsequent iterations we could accept windowStart+1, windowStart+2, etc.
                // So the only region we can confidently reject is this current position.
                rejectedPositions.add(windowStart);
                ++windowStart;
            }
        }

        // Compute merged rejected regions.
        List<RejectedRegion> rejectedRegions = new ArrayList<>();
        if(!rejectedPositions.isEmpty())
        {
            String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
            int rejectStart = rejectedPositions.get(0);
            int rejectEnd = rejectedPositions.get(0);
            for(int i = 1; i < rejectedPositions.size(); ++i)
            {
                int rejectedPosition = rejectedPositions.get(i);
                if(rejectedPosition != rejectEnd + 1 || i == rejectedPositions.size() - 1)
                {
                    rejectedRegions.add(
                            new RejectedRegion(new ChrBaseRegion(target.chromosome(), rejectStart, rejectEnd), source, rejectionReason));
                    rejectStart = rejectedPosition;
                }
                rejectEnd = rejectedPosition;
            }
        }

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
        return outwardMovingCenteredProbes(new BasePosition(target.chromosome(), regionCentre), target, source)
                // Stop once the probes are too far from the target position to cover it.
                .takeWhile(probe -> probeCoversTarget(probe.probeRegion().baseRegion(), target.baseRegion()));
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
                .takeWhile(probe -> probe.probeRegion().containsPosition(position.Position));
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
        return new CandidateProbe(source, new ChrBaseRegion(chromosome, startPosition, startPosition + PROBE_LENGTH - 1), targetRegion);
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
            case MAX_QUALITY -> selectBestScoringProbe(acceptableProbes,
                    EvaluatedProbe::qualityScore,
                    quality -> Doubles.greaterOrEqual(quality, 1.0), true);
            // Stop if a probe with gc=gcTarget is found.
            case BEST_GC -> selectBestScoringProbe(acceptableProbes,
                    probe -> abs(probe.gcContent() - criteria.eval().gcContentTarget()),
                    gc -> Doubles.equal(gc, criteria.eval().gcContentTarget()), false);
        };
    }

    // Min/max function with early stopping if an optimal value is found.
    private static Optional<EvaluatedProbe> selectBestScoringProbe(Stream<EvaluatedProbe> probes,
            final ToDoubleFunction<EvaluatedProbe> scoreFunc, final DoublePredicate isOptimalFunc, boolean maximise)
    {
        BiPredicate<Double, Double> scoreCompareFunc = maximise ? Doubles::greaterThan : Doubles::lessThan;
        Optional<EvaluatedProbe> bestProbe = Optional.empty();
        double bestScore = maximise ? NEGATIVE_INFINITY : POSITIVE_INFINITY;
        Iterator<EvaluatedProbe> iterator = probes.iterator();
        while(iterator.hasNext())
        {
            EvaluatedProbe probe = iterator.next();
            double score = scoreFunc.applyAsDouble(probe);
            if(scoreCompareFunc.test(score, bestScore))
            {
                bestProbe = Optional.of(probe);
                bestScore = score;
                LOGGER.trace("New bestScore={} probe={}", bestScore, bestProbe);
                if(isOptimalFunc.test(bestScore))
                {
                    LOGGER.trace("Optimal score found, stopping candidate search");
                    break;
                }
            }
        }
        return bestProbe;
    }
}
