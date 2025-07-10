package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.function.Consumer;

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
    // Hook to catch all candidate probes for output.
    private final Consumer<EvaluatedProbe> mCandidateCallback;

    private static final Logger LOGGER = LogManager.getLogger(ProbeEvaluator.class);

    public ProbeEvaluator(final RefGenomeInterface refGenome, final ProbeQualityProfile qualityProfile,
            final Consumer<EvaluatedProbe> candidateCallback)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
        mCandidateCallback = candidateCallback;
    }

    public EvaluatedProbe evaluateCandidate(final CandidateProbe probe, final ProbeEvalCriteria criteria)
    {
        EvaluatedProbe evaluatedProbe = new EvaluatedProbe(probe, criteria);
        evaluatedProbe = evaluateQualityScore(evaluatedProbe, criteria);
        if(!evaluatedProbe.rejected())
        {
            evaluatedProbe = evaluateGcContent(evaluatedProbe, criteria);
        }
        logCandidateProbe(evaluatedProbe);
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
            // Maybe be interesting to know when this happens because the probe quality profile ideally covers the whole genome.
            LOGGER.trace("Candidate probe not covered by probe quality profile so assuming qualityScore=0 probe={}", probe);
            return 0d;
        });
    }

    private String getProbeSequence(final CandidateProbe probe)
    {
        ChrBaseRegion region = probe.probeRegion();
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }

    private void logCandidateProbe(final EvaluatedProbe probe)
    {
        LOGGER.trace("Evaluated probe: {}", probe);
        mCandidateCallback.accept(probe);
    }
}
