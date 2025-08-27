package com.hartwig.hmftools.cobalt.utils;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public record RawCobaltRatio(
        @NotNull String chromosome,
        int position,
        double referenceReadCount,
        double tumorReadCount,
        double referenceGcRatio,
        double tumorGcRatio,
        double referenceGcDiploidRatio
)
{
    public boolean matchesPosition(final RawCobaltRatio other)
    {
        return chromosome.equals(other.chromosome) && position == other.position;
    }

    public RawCobaltRatio differences(final RawCobaltRatio other, final double epsilon)
    {
        List<DoubleDifference> differences = Lists.newArrayList();
        differences.add(new DoubleDifference(referenceReadCount, other.referenceReadCount, epsilon));
        differences.add(new DoubleDifference(tumorReadCount, other.tumorReadCount, epsilon));
        differences.add(new DoubleDifference(referenceGcRatio, other.referenceGcRatio, epsilon));
        differences.add(new DoubleDifference(tumorGcRatio, other.tumorGcRatio, epsilon));
        differences.add(new DoubleDifference(referenceGcDiploidRatio, other.referenceGcDiploidRatio, epsilon));
        boolean hasDifference = differences.stream().anyMatch(diff -> diff.hasDifference);
        if(hasDifference)
        {
            return new RawCobaltRatio(chromosome, position, differences.get(0).difference, differences.get(1).difference, differences.get(2).difference, differences.get(3).difference, differences.get(4).difference);
        }
        return null;
    }
}
