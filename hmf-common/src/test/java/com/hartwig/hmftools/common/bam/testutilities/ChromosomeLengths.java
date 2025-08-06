package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public interface ChromosomeLengths
{
    int chromosomeLength(String chromosomeName);
}

class RefGenomeBackedChromosomeLengths implements ChromosomeLengths
{
    private final RefGenomeInterface refGenomeInterface;

    RefGenomeBackedChromosomeLengths(final RefGenomeInterface refGenomeInterface)
    {
        this.refGenomeInterface = refGenomeInterface;
    }

    @Override
    public int chromosomeLength(final String chromosomeName)
    {
        return refGenomeInterface.getChromosomeLength(chromosomeName);
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
    public int chromosomeLength(final String chromosomeName)
    {
        return length;
    }
}
