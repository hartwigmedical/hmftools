package com.hartwig.hmftools.common.bam.testutilities;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public interface ChromosomeLengths
{
    int chromosomeLength(HumanChromosome chromosomeName);
}

class RefGenomeBackedChromosomeLengths implements ChromosomeLengths
{
    private final RefGenomeInterface refGenomeInterface;

    RefGenomeBackedChromosomeLengths(final RefGenomeInterface refGenomeInterface)
    {
        this.refGenomeInterface = refGenomeInterface;
    }

    @Override
    public int chromosomeLength(final HumanChromosome chromosome)
    {
        return refGenomeInterface.getChromosomeLength(V38.versionedChromosome(chromosome));
    }
}

class ConstantChromosomeLengths implements ChromosomeLengths
{
    private final int length;

    ConstantChromosomeLengths(final int length)
    {
        this.length = length;
    }

    @Override
    public int chromosomeLength(final HumanChromosome chromosome)
    {
        return length;
    }
}
