package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;
import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.DoubleSupplier;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
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

    // Creates a probe from its target region(s)/sequence.
    // Returns empty optional if it's not valid to create a probe at that location.
    public Optional<Probe> createProbe(final ProbeTarget target, final TargetMetadata metadata)
    {
        String sequence = buildSequence(target);
        if(sequence.length() != target.baseLength())
        {
            return Optional.empty();
        }

        if(target.isExactRegion())
        {
            return createProbe(
                    target, sequence, metadata,
                    () -> getQualityScore(target.exactRegion(), sequence, false),
                    () -> calcGcPercent(sequence));
        }
        else
        {
            return createProbe(
                    target, sequence, metadata,
                    () -> getQualityScore(null, sequence, true),
                    () -> calcGcPercent(sequence));
        }
    }

    private Optional<Probe> createProbe(final ProbeTarget target, final String sequence, final TargetMetadata metadata,
            final DoubleSupplier getQualityScore, final DoubleSupplier getGcContent)
    {
        // Only check properties which are inconvenient for the caller to check in advance.
        // Everything else is expected to be checked by the caller and will generate an exception.
        boolean regionsValid = target.regions().stream().allMatch(this::isRegionValid);
        boolean sequenceValid = isDnaSequenceNormal(sequence);
        boolean valid = regionsValid && sequenceValid;
        if(valid)
        {
            return Optional.of(new Probe(
                    target, sequence, metadata,
                    null, null,
                    getQualityScore.getAsDouble(), getGcContent.getAsDouble()));
        }
        else
        {
            return Optional.empty();
        }
    }

    private boolean isRegionValid(final ChrBaseRegion region)
    {
        return region.hasValidPositions() && region.end() <= mRefGenome.getChromosomeLength(region.chromosome());
    }

    private String buildSequence(final ProbeTarget target)
    {
        String start = target.startRegion() == null ? "" : getSequence(target.startRegion());
        if(target.startOrientation() == Orientation.REVERSE)
        {
            start = reverseComplementBases(start);
        }
        String insert = target.insertSequence() == null ? "" : target.insertSequence();
        String end = target.endRegion() == null ? "" : getSequence(target.endRegion());
        if(target.endOrientation() == Orientation.REVERSE)
        {
            end = reverseComplementBases(end);
        }
        return start + insert + end;
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
