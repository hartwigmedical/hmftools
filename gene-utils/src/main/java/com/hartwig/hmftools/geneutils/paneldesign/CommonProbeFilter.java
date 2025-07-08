package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Probe candidates filtering that applies to probes from all sources.
public class CommonProbeFilter
{
    private final RefGenomeInterface mRefGenome;
    private final ProbeQualityProfile mQualityProfile;
    private final double mQualityMin;
    private final double mGcMin;
    private final double mGcMax;

    public CommonProbeFilter(RefGenomeInterface refGenome, final ProbeQualityProfile qualityProfile,
            double qualityMin, double gcMin, double gcMax)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
        mQualityMin = qualityMin;
        mGcMin = gcMin;
        mGcMax = gcMax;
    }

    // Returns filter reason if filter fails, otherwise empty optional.
    public Optional<String> filter(final ProbeCandidate probe)
    {
        Stream<Function<ProbeCandidate, Optional<String>>> filters = Stream.of(
                this::filterQuality, this::filterGc);
        return filters
                .map(filter -> filter(probe))
                .filter(Optional::isPresent)
                .findFirst()
                .orElse(Optional.empty());
    }

    private Optional<String> filterQuality(final ProbeCandidate probe)
    {
        Optional<Double> qualityScore = mQualityProfile.computeQualityScore(probe.probeRegion());
        if(qualityScore.map(q -> q >= mQualityMin).orElse(false))
        {
            return Optional.empty();
        }
        else
        {
            return Optional.of("quality score");
        }
    }

    private Optional<String> filterGc(final ProbeCandidate probe)
    {
        ChrBaseRegion region = probe.probeRegion();
        String sequence = mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
        double gcContent = calcGcPercent(sequence);
        if(gcContent >= mGcMin && gcContent <= mGcMax)
        {
            return Optional.empty();
        }
        else
        {
            return Optional.of("gc");
        }
    }
}
