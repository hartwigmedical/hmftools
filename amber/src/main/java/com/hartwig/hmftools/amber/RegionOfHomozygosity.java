package com.hartwig.hmftools.amber;

import java.util.Comparator;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class RegionOfHomozygosity implements Comparable<RegionOfHomozygosity>
{
    public final Chromosome Chromosome;
    public final int Start;
    public final int End;
    public final int NumHomozygous;
    public final int NumHeterozygous;
    public final int NumUnclear;

    public int getSnpCount() { return NumHomozygous + NumHeterozygous + NumUnclear; }
    public int getLength() { return End - Start + 1; }

    public RegionOfHomozygosity(final Chromosome chromosome, int start, int end, int numHomozygous, int numHeterozygous, int numUnclear)
    {
        Chromosome = chromosome;
        Start = start;
        End = end;
        NumHomozygous = numHomozygous;
        NumHeterozygous = numHeterozygous;
        NumUnclear = numUnclear;
    }

    @Override
    public int compareTo(final RegionOfHomozygosity other)
    {
        int chrRank = HumanChromosome.chromosomeRank(Chromosome.toString());
        int chrRankOther = HumanChromosome.chromosomeRank(other.Chromosome.toString());

        if(chrRank != chrRankOther)
            return chrRank < chrRankOther ? -1 : 1;

        if(Start != other.Start)
            return Start < other.Start ? -1 : 1;

        if(End != other.End)
            return End < other.End ? -1 : 1;

        return 0;
    }

}
