package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.utils.Doubles;

public record CobaltChromosome(HumanChromosome humanChromosome, double typicalRatio, double actualRatio, boolean mosiac)
        implements Chromosome
{
    public double typicalRatio()
    {
        return typicalRatio;
    }

    public double actualRatio()
    {
        return actualRatio;
    }

    public boolean mosiac()
    {
        return mosiac;
    }

    public boolean isNormal()
    {
        return Doubles.equal(typicalRatio(), actualRatio());
    }

    @Override
    public boolean isAutosome()
    {
        return humanChromosome().isAutosome();
    }

    @Override
    public boolean isAllosome()
    {
        return humanChromosome().isAllosome();
    }

    public boolean isDiploid()
    {
        return Doubles.equal(actualRatio(), 1.0);
    }
}
