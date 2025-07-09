package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.Comparator;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

    // TODO: if no best probe, return reason why
    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there is no such best probe.
    public Optional<EvaluatedProbe> selectBestProbe(Stream<CandidateProbe> probes, final ProbeSelectCriteria criteria)
    {
        Stream<EvaluatedProbe> acceptableProbes = probes
                .map(probe -> evaluateCandidate(probe, criteria.eval()))
                .filter(EvaluatedProbe::accepted);
        return switch(criteria.strategy())
        {
            case FIRST_ACCEPTABLE -> acceptableProbes.findFirst();
            case MAX_QUALITY -> acceptableProbes.max(Comparator.nullsFirst(Comparator.comparing(EvaluatedProbe::qualityScore)));
            case BEST_GC -> acceptableProbes
                    .min(Comparator.comparingDouble(
                            probe -> abs(probe.gcContent() - criteria.eval().gcContentTarget())));
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
}
