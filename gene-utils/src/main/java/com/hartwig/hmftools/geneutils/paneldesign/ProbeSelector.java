package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.geneutils.paneldesign.Utils.getBestScoringElement;

import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.Doubles;

public class ProbeSelector
{
    private final ProbeEvaluator mProbeEvaluator;

    public ProbeSelector(final ProbeEvaluator probeEvaluator)
    {
        mProbeEvaluator = probeEvaluator;
    }

    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there are no acceptable probes.
    public Optional<EvaluatedProbe> selectBestProbe(Stream<CandidateProbe> probes, final ProbeSelectCriteria criteria)
    {
        Stream<EvaluatedProbe> acceptableProbes = probes
                .map(probe -> mProbeEvaluator.evaluateCandidate(probe, criteria.eval()))
                .filter(EvaluatedProbe::accepted);

        ProbeSelectStrategy strategy = criteria.strategy();
        if(strategy instanceof ProbeSelectStrategy.FirstAcceptable)
        {
            return acceptableProbes.findFirst();
        }
        else if(strategy instanceof ProbeSelectStrategy.MaxQuality)
        {
            // Early stopping if "optimal" quality score is found.
            double optimalQuality = ((ProbeSelectStrategy.MaxQuality) strategy).optimalQuality();
            return getBestScoringElement(
                    acceptableProbes,
                    EvaluatedProbe::qualityScore,
                    quality -> Doubles.greaterOrEqual(quality, optimalQuality),
                    true);
        }
        else if(strategy instanceof ProbeSelectStrategy.BestGc)
        {
            // Early stopping if "optimal" GC content is found.
            double targetGc = criteria.eval().gcContentTarget();
            double optimalGcTolerance = ((ProbeSelectStrategy.BestGc) strategy).gcToleranceOptimal();
            return getBestScoringElement(
                    acceptableProbes,
                    probe -> abs(probe.gcContent() - targetGc),
                    distance -> Doubles.lessOrEqual(distance, optimalGcTolerance),
                    false);
        }
        else
        {
            throw new IllegalArgumentException("Unhandled ProbeSelectStrategy");
        }
    }
}
