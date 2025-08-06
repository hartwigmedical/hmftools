package com.hartwig.hmftools.cobalt.e2e;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public record ConstantGcFileSection(HumanChromosome chromosome, int start, int stop, double gcRatio) implements GcFileSection
{
    @Override
    public double gcRatio(final int position)
    {
        return gcRatio;
    }
}
