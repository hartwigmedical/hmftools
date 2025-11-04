package com.hartwig.hmftools.common.cobalt;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public record CobaltRatio(
        @NotNull String chromosome,
        int position,
        double referenceReadDepth,
        double referenceGCRatio,
        double referenceGcContent,
        double referenceGCDiploidRatio,
        double tumorReadDepth,
        double tumorGCRatio,
        double tumorGcContent
) implements GenomePosition
{
    public static final double COBALT_MASKED_VALUE = -1.0;

    public CobaltRatio realign(int newPosition)
    {
        return new CobaltRatio(chromosome, newPosition, referenceReadDepth, referenceGCRatio, referenceGcContent,
                referenceGCDiploidRatio, tumorReadDepth, tumorGCRatio, tumorGcContent);
    }

    public CobaltRatio normaliseTumorGcRatio(double factor)
    {
        return new CobaltRatio(chromosome, position, referenceReadDepth, referenceGCRatio, referenceGcContent,
                referenceGCDiploidRatio, tumorReadDepth, tumorGCRatio / factor, tumorGcContent);
    }

    public CobaltRatio mask()
    {
        return new CobaltRatio(
                chromosome, position, referenceReadDepth, COBALT_MASKED_VALUE, COBALT_MASKED_VALUE, COBALT_MASKED_VALUE,
                tumorReadDepth, COBALT_MASKED_VALUE, COBALT_MASKED_VALUE);
    }

    public ChrBaseRegion window()
    {
        return new ChrBaseRegion(chromosome, position, position + 1000 - 1);
    }

    public <T extends ChrBaseRegion> List<T> findWindowOverlaps(@NotNull List<T> intervals)
    {
        return intervals.stream().filter(t -> this.window().overlaps(t)).collect(Collectors.toList());
    }

    @Override
    public String toString()
    {
        return "CobaltRatio(" + chromosome + ", " + position + ")";
    }
}
