package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;

import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

import org.jetbrains.annotations.Nullable;

// Helps compute probe quality scores.
// Takes quality scores from the probe quality profile, falling back to the alignment model.
// Also batches computations on the alignment model to make it faster.
public class ProbeQualityScorer
{
    private final ProbeQualityProfile mQualityProfile;
    private final ProbeQualityModel mQualityModel;

    public ProbeQualityScorer(final ProbeQualityProfile qualityProfile, final ProbeQualityModel qualityModel)
    {
        mQualityProfile = qualityProfile;
        mQualityModel = qualityModel;
    }

    public List<Probe> getQualityScore(Stream<Probe> probes)
    {
        // TODO: batch
        return probes
                .map(probe -> probe.withQualityScore(
                        getQualityScore(probe.definition().exactRegionOrNull(), probe.sequence(), !probe.definition().isExactRegion())))
                .toList();
    }

    private double getQualityScore(@Nullable final ChrBaseRegion region, @Nullable final String sequence, boolean allowFallback)
    {
        // Prefer the quality profile first since it's much faster.
        // Fall back to computing the quality score based on the sequence if required.

        OptionalDouble qualityScore = mQualityProfile != null && region != null
                ? mQualityProfile.computeQualityScore(region)
                : OptionalDouble.empty();
        if(qualityScore.isEmpty() && mQualityModel != null && sequence != null && allowFallback)
        {
            ProbeQualityModel.Result modelResult = mQualityModel.compute(List.of(sequence.getBytes())).get(0);
            qualityScore = OptionalDouble.of(modelResult.qualityScore());
        }
        return qualityScore.orElse(DEFAULT_PROBE_QUALITY);
    }
}
