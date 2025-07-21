package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.geneutils.paneldesign.Utils.getBestScoringElement;

import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.Doubles;

// Evaluates and selects best candidate probes.
// TODO: review interface. should this contain ProbeEvaluator?
public class ProbeSelector
{
    public final ProbeEvaluator mProbeEvaluator;

    public ProbeSelector(final ProbeEvaluator probeEvaluator)
    {
        mProbeEvaluator = probeEvaluator;
    }

    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there are no acceptable probes.
    public Optional<EvaluatedProbe> selectBestCandidate(Stream<CandidateProbe> probes, final Criteria criteria)
    {
        Stream<EvaluatedProbe> evaluatedProbes = probes.map(probe -> mProbeEvaluator.evaluateCandidate(probe, criteria.eval()));
        return selectBestProbe(evaluatedProbes, criteria);
    }

    private Optional<EvaluatedProbe> selectBestProbe(Stream<EvaluatedProbe> probes, final Criteria criteria)
    {
        Stream<EvaluatedProbe> acceptableProbes = probes.filter(EvaluatedProbe::accepted);

        Strategy strategy = criteria.strategy();
        if(strategy instanceof Strategy.FirstAcceptable)
        {
            return acceptableProbes.findFirst();
        }
        else if(strategy instanceof Strategy.MaxQuality)
        {
            // Early stopping if "optimal" quality score is found.
            double optimalQuality = ((Strategy.MaxQuality) strategy).optimalQuality();
            return getBestScoringElement(
                    acceptableProbes,
                    EvaluatedProbe::qualityScore,
                    quality -> Doubles.greaterOrEqual(quality, optimalQuality),
                    true);
        }
        else if(strategy instanceof Strategy.BestGc)
        {
            // Early stopping if "optimal" GC content is found.
            double targetGc = criteria.eval().gcContentTarget();
            double optimalGcTolerance = ((Strategy.BestGc) strategy).gcToleranceOptimal();
            return getBestScoringElement(
                    acceptableProbes,
                    probe -> abs(probe.gcContent() - targetGc),
                    distance -> Doubles.lessOrEqual(distance, optimalGcTolerance),
                    false);
        }
        else
        {
            throw new IllegalArgumentException("Unhandled Strategy");
        }
    }

    // When there are multiple acceptable candidate probes, how to select the best?
    public sealed interface Strategy
            permits
            Strategy.FirstAcceptable, Strategy.MaxQuality, Strategy.BestGc
    {
        // Pick the first probe that is acceptable.
        record FirstAcceptable() implements Strategy
        {
        }

        // Pick the acceptable probe with the highest quality score.
        record MaxQuality(
                // Consider the max quality to have been found if exceeding this value.
                // Useful to improve runtime performance.
                double optimalQuality
        ) implements Strategy
        {
            public MaxQuality()
            {
                this(1);
            }
        }

        // Pick the acceptable probe with the GC content closest to gcContentTarget.
        record BestGc(
                // Consider the best GC to have been found if within this tolerance.
                // Useful to improve runtime performance.
                double gcToleranceOptimal
        ) implements Strategy
        {
        }
    }

    public record Criteria(
            ProbeEvaluator.Criteria eval,
            Strategy strategy
    )
    {
    }
}
