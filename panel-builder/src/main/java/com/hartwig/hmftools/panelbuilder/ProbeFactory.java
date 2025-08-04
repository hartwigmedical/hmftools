package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;

import java.util.List;
import java.util.OptionalDouble;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class ProbeFactory
{
    private final RefGenomeInterface mRefGenome;
    @Nullable
    private final ProbeQualityProfile mQualityProfile;

    private static final Logger LOGGER = LogManager.getLogger(ProbeFactory.class);

    public ProbeFactory(final RefGenomeInterface refGenome, @Nullable final ProbeQualityProfile qualityProfile)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
    }

    // Create a probe that corresponds exactly to a region in the reference genome.
    public Probe createProbeFromRegion(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        String sequence = getSequence(region);
        return new Probe(
                region, sequence, metadata,
                null, null,
                getQualityScore(region), calcGcPercent(sequence));
    }

    // Create a probe with a custom sequence that's built from multiple regions in the reference genome.
    public Probe createProbeFromSequence(final String sequence, final TargetMetadata metadata,
            final List<ChrBaseRegion> contributingRegions)
    {
        double qualityScore = contributingRegions.stream().mapToDouble(this::getQualityScore).min().orElseThrow();
        return new Probe(null, sequence, metadata, null, null, qualityScore, calcGcPercent(sequence));
    }

    private String getSequence(final ChrBaseRegion region)
    {
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }

    private double getQualityScore(final ChrBaseRegion region)
    {
        // Allow mQualityProfile to be null for testing purposes where the quality score doesn't matter.
        OptionalDouble qualityScore = mQualityProfile == null ? OptionalDouble.empty() : mQualityProfile.computeQualityScore(region);
        return qualityScore.orElseGet(() ->
        {
            double quality = DEFAULT_PROBE_QUALITY;
            LOGGER.trace("Probe region not covered by probe quality profile so assuming qualityScore={} probe={}",
                    quality, region);
            return quality;
        });
    }
}
