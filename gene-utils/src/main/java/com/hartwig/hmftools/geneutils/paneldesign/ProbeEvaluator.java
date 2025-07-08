package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Common candidate probe evaluation and filtering.
public class ProbeEvaluator
{
    private final RefGenomeInterface mRefGenome;
    private final ProbeQualityProfile mQualityProfile;
    private final double mQualityMin;
    private final double mGcMin;
    private final double mGcMax;

    public ProbeEvaluator(RefGenomeInterface refGenome, final ProbeQualityProfile qualityProfile,
            double qualityMin, double gcMin, double gcMax)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
        mQualityMin = qualityMin;
        mGcMin = gcMin;
        mGcMax = gcMax;
    }

    // Gets the best acceptable probe from a set of candidate probes. Returns empty optional if there is no such best probe.
    public Optional<EvaluatedProbe> selectBestProbe(Stream<CandidateProbe> probes)
    {
        // Apply the common probe filters.
        // Select only probes which pass the filters.
        // Sort by highest quality first.
        // Take the first (i.e. best) acceptable probe.
        return probes
                // Apply the common probe filters.
                .map(this::evaluateCandidate)
                // Select only probes which pass the filters.
                .filter(EvaluatedProbe::accepted)
                // Order nulls as lower values to be safe, but a quality score is required so it shouldn't matter here.
                .max(Comparator.nullsFirst(Comparator.comparing(EvaluatedProbe::qualityScore)));
    }

    public EvaluatedProbe evaluateCandidate(final CandidateProbe probe)
    {
        List<Function<EvaluatedProbe, EvaluatedProbe>> filters = List.of(this::evaluateQuality, this::evaluateGc);
        EvaluatedProbe evaluatedProbe = new EvaluatedProbe(probe);
        for(Function<EvaluatedProbe, EvaluatedProbe> filter : filters)
        {
            evaluatedProbe = filter.apply(evaluatedProbe);
            if(evaluatedProbe.rejected())
            {
                // No need to keep going if a filter already rejected the probe.
                break;
            }
        }
        return evaluatedProbe;
    }

    private EvaluatedProbe evaluateQuality(EvaluatedProbe probe)
    {
        double qualityScore = getProbeQuality(probe.candidate());
        probe = probe.withQualityScore(qualityScore);
        if(!(qualityScore >= mQualityMin))
        {
            probe = probe.withRejectionReason("quality");
        }
        return probe;
    }

    private EvaluatedProbe evaluateGc(EvaluatedProbe probe)
    {
        double gcContent = getGcContent(probe.candidate());
        probe = probe.withGcContent(gcContent);
        if(!(gcContent >= mGcMin && gcContent <= mGcMax))
        {
            probe = probe.withRejectionReason("gc");
        }
        return probe;
    }

    private double getProbeQuality(CandidateProbe probe)
    {
        // Never want to accept a probe with no quality score, so just return 0 in that case to simplify the code elsewhere.
        return mQualityProfile.computeQualityScore(probe.probeRegion()).orElse(0d);
    }

    private double getGcContent(CandidateProbe probe)
    {
        ChrBaseRegion region = probe.probeRegion();
        String sequence = mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
        double gcContent = calcGcPercent(sequence);
        return gcContent;
    }
}
