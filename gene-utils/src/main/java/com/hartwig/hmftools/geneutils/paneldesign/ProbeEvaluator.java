package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Common probe candidate evaluation and filtering.
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
    public Optional<ProbeCandidate> selectBestProbe(Stream<ProbeCandidate> probes)
    {
        // Apply the common probe filters.
        // Select only probes which pass the filters.
        // Sort by highest quality first.
        // Take the first (i.e. best) acceptable probe.
        return probes
                // Apply the common probe filters.
                .map(this::filter)
                // Select only probes which pass the filters.
                .filter(EvaluatedProbeCandidate::acceptable)
                .min(Comparator.nullsLast(Comparator.comparing(probe -> probe.QualityScore)))
                .map(probe -> probe.Probe);
    }

    public static class EvaluatedProbeCandidate
    {
        public final ProbeCandidate Probe;
        // null if probe is acceptable.
        public String RejectionReason = null;
        // null if not checked.
        public Double QualityScore = null;
        // null if not checked.
        public Double GcContent = null;

        EvaluatedProbeCandidate(final ProbeCandidate probe)
        {
            Probe = probe;
        }

        public boolean acceptable()
        {
            return RejectionReason == null;
        }
    }

    public EvaluatedProbeCandidate filter(final ProbeCandidate probe)
    {
        List<Consumer<EvaluatedProbeCandidate>> filters = List.of(this::filterQuality, this::filterGc);

        EvaluatedProbeCandidate filteredProbe = new EvaluatedProbeCandidate(probe);
        for(Consumer<EvaluatedProbeCandidate> filter : filters)
        {
            filter.accept(filteredProbe);
            if(!filteredProbe.acceptable())
            {
                break;
            }
        }
        return filteredProbe;
    }

    private void filterQuality(EvaluatedProbeCandidate probe)
    {
        double qualityScore = getProbeQuality(probe.Probe);
        probe.QualityScore = qualityScore;
        if(!(qualityScore >= mQualityMin))
        {
            probe.RejectionReason = "quality";
        }
    }

    private void filterGc(EvaluatedProbeCandidate probe)
    {
        double gcContent = getGcContent(probe.Probe);
        probe.GcContent = gcContent;
        if(!(gcContent >= mGcMin && gcContent <= mGcMax))
        {
            probe.RejectionReason = "gc";
        }
    }

    private double getProbeQuality(ProbeCandidate probe)
    {
        // Never want to accept a probe with no quality score, so just return 0 in that case to simplify the code elsewhere.
        return mQualityProfile.computeQualityScore(probe.probeRegion()).orElse(0d);
    }

    private double getGcContent(ProbeCandidate probe)
    {
        ChrBaseRegion region = probe.probeRegion();
        String sequence = mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
        double gcContent = calcGcPercent(sequence);
        return gcContent;
    }
}
