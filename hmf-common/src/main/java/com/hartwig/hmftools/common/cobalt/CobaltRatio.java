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

    public CobaltRatio realign(int newPosition)
    {
        return new CobaltRatio(chromosome, newPosition, referenceReadDepth, referenceGCRatio, referenceGcContent,
                referenceGCDiploidRatio, tumorReadDepth, tumorGCRatio, tumorGcContent);
    }

    public CobaltRatio mask()
    {
        return new CobaltRatio(chromosome, position, referenceReadDepth, -1.0, -1.0, -1.0, tumorReadDepth, -1.0, -1.0);
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
