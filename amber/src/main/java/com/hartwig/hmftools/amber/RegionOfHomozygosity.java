package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class RegionOfHomozygosity
{
    private final Chromosome mChromosome;
    private final int mStart;
    private final int mEnd;
    private final int mNumHomozygous;
    private final int mNumHeterozygous;
    private final int mNumUnclear;

    public Chromosome getChromosome() { return mChromosome; }
    public int getStart() { return mStart; }
    public int getEnd() { return mEnd; }
    public int getNumHomozygous() { return mNumHomozygous; }
    public int getNumHeterozygous() { return mNumHeterozygous; }
    public int getNumUnclear() { return mNumUnclear; }
    public int getSnpCount() { return mNumHomozygous + mNumHeterozygous + mNumUnclear; }
    public int getLength() { return mEnd - mStart + 1; }

    public RegionOfHomozygosity(Chromosome chromosome, int start, int end, int numHomozygous, int numHeterozygous, int numUnclear)
    {
        mChromosome = chromosome;
        mStart = start;
        mEnd = end;
        mNumHomozygous = numHomozygous;
        mNumHeterozygous = numHeterozygous;
        mNumUnclear = numUnclear;
    }
}
