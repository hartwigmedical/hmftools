package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.Iterator;
import java.util.Optional;
import java.util.function.BiPredicate;
import java.util.function.DoublePredicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Common candidate probe evaluation and filtering.
public class ProbeEvaluator
{
    private final RefGenomeInterface mRefGenome;
    private final ProbeQualityProfile mQualityProfile;

    private static final Logger LOGGER = LogManager.getLogger(ProbeEvaluator.class);

    public ProbeEvaluator(final RefGenomeInterface refGenome, final ProbeQualityProfile qualityProfile)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
    }

    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there are no acceptable probes.
    public Optional<EvaluatedProbe> selectBestProbe(Stream<CandidateProbe> probes, final ProbeSelectCriteria criteria)
    {
        Stream<EvaluatedProbe> acceptableProbes = probes
                .map(probe -> evaluateCandidate(probe, criteria.eval()))
                .filter(EvaluatedProbe::accepted);
        return switch(criteria.strategy())
        {
            case FIRST_ACCEPTABLE -> acceptableProbes.findFirst();
            // Stop if a probe with quality=1 is found.
            case MAX_QUALITY -> findBestProbe(acceptableProbes,
                    EvaluatedProbe::qualityScore,
                    quality -> Doubles.greaterOrEqual(quality, 1.0), true);
            // Stop if a probe with gc=gcTarget is found.
            case BEST_GC -> findBestProbe(acceptableProbes,
                    probe -> abs(probe.gcContent() - criteria.eval().gcContentTarget()),
                    gc -> Doubles.equal(gc, criteria.eval().gcContentTarget()), false);
        };
    }

    private EvaluatedProbe evaluateCandidate(final CandidateProbe probe, final ProbeEvalCriteria criteria)
    {
        EvaluatedProbe evaluatedProbe = new EvaluatedProbe(probe, criteria);
        evaluatedProbe = evaluateQualityScore(evaluatedProbe, criteria);
        if(!evaluatedProbe.rejected())
        {
            evaluatedProbe = evaluateGcContent(evaluatedProbe, criteria);
        }
        LOGGER.trace("Evaluated probe: {}", evaluatedProbe);
        return evaluatedProbe;
    }

    private EvaluatedProbe evaluateQualityScore(EvaluatedProbe probe, final ProbeEvalCriteria criteria)
    {
        double qualityScore = getProbeQuality(probe.candidate());
        probe = probe.withQualityScore(qualityScore);
        if(!(qualityScore >= criteria.qualityScoreMin()))
        {
            probe = probe.withRejectionReason("quality");
        }
        return probe;
    }

    private EvaluatedProbe evaluateGcContent(EvaluatedProbe probe, final ProbeEvalCriteria criteria)
    {
        String sequence = getProbeSequence(probe.candidate());
        probe = probe.withSequence(sequence);
        double gcContent = calcGcPercent(sequence);
        probe = probe.withGcContent(gcContent);
        if(!(gcContent >= criteria.gcContentMin() && gcContent <= criteria.gcContentMax()))
        {
            probe = probe.withRejectionReason("gc");
        }
        return probe;
    }

    private double getProbeQuality(final CandidateProbe probe)
    {
        return mQualityProfile.computeQualityScore(probe.probeRegion()).orElseGet(() ->
        {
            // Never want to accept a probe with no quality score, so just return 0 in that case to simplify the code elsewhere.
            // But also this shouldn't happen in practice because the probe quality profile covers the whole genome.
            LOGGER.warn("Evaluating candidate probe not covered by probe quality profile; assuming qualityScore=0");
            return 0d;
        });
    }

    private String getProbeSequence(final CandidateProbe probe)
    {
        ChrBaseRegion region = probe.probeRegion();
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }

    // Min/max function with early stopping if an optimal value is found.
    private static Optional<EvaluatedProbe> findBestProbe(Stream<EvaluatedProbe> probes, final ToDoubleFunction<EvaluatedProbe> scoreFunc,
            final DoublePredicate isOptimalFunc, boolean maximise)
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
