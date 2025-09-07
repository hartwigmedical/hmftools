package com.hartwig.hmftools.cobalt.e2e;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public record ConstantMappablePercentageGcFileSection(HumanChromosome chromosome, int start, int stop, double mappablePercentage) implements GcFileSection
{
    @Override
    public double gcRatio(final int position)
    {
        return 0.5;
    }

    @Override
    public double mappablePercentage(final int position)
    {
        return mappablePercentage;
    }
}
