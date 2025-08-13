package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;
import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.DoubleSupplier;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

import org.jetbrains.annotations.Nullable;

// Handles creation of individual probes.
public class ProbeFactory
{
    private final RefGenomeInterface mRefGenome;
    // Allow these to be null for testing purposes where the quality score doesn't matter.
    @Nullable
    private final ProbeQualityProfile mQualityProfile;
    @Nullable
    private final ProbeQualityModel mQualityModel;

    public ProbeFactory(final RefGenomeInterface refGenome, @Nullable final ProbeQualityProfile qualityProfile,
            @Nullable final ProbeQualityModel qualityModel)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
        mQualityModel = qualityModel;
    }

    // Create a probe that corresponds exactly to a region in the reference genome.
    // Returns empty optional if it's not valid to create a probe at that location.
    public Optional<Probe> createProbeFromRegion(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        String sequence = getSequence(region);
        return createProbe(
                region, sequence, metadata,
                () -> getQualityScore(region, sequence, false),
                () -> calcGcPercent(sequence));
    }

    // Create a probe with a custom sequence.
    // Returns empty optional if it's not valid to create a probe at that location.
    public Optional<Probe> createProbeFromSequence(final String sequence, final TargetMetadata metadata)
    {
        return createProbe(
                null, sequence, metadata,
                () -> getQualityScore(null, sequence, true),
                () -> calcGcPercent(sequence));
    }

    private Optional<Probe> createProbe(@Nullable final ChrBaseRegion region, final String sequence, final TargetMetadata metadata,
            final DoubleSupplier getQualityScore, final DoubleSupplier getGcContent)
    {
        // Only check properties which are inconvenient for the caller to check in advance.
        // Everything else is expected to be checked by the caller and will generate an exception.
        boolean regionValid =
                region == null || (region.start() >= 1 && region.end() <= mRefGenome.getChromosomeLength(region.chromosome()));
        boolean sequenceValid = isDnaSequenceNormal(sequence);
        boolean valid = regionValid && sequenceValid;
        if(valid)
        {
            return Optional.of(new Probe(
                    region, sequence, metadata,
                    null, null,
                    getQualityScore.getAsDouble(), getGcContent.getAsDouble()));
        }
        else
        {
            return Optional.empty();
        }
    }

    private String getSequence(final ChrBaseRegion region)
    {
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
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
