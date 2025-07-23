package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

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

    public CobaltRatio realign(int newPosition)
    {
        return new CobaltRatio(chromosome, newPosition, referenceReadDepth, referenceGCRatio, referenceGcContent,
                referenceGCDiploidRatio, tumorReadDepth, tumorGCRatio, tumorGcContent);
    }
}
