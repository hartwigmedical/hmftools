package com.hartwig.hmftools.common.utils.sv;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public final class SvCommonUtils
{
    public static final byte POS_ORIENT = 1;
    public static final byte NEG_ORIENT = -1;

    public static boolean lowerChromosome(final String chr, final String otherChr)
    {
        return chromosomeRank(chr) < chromosomeRank(otherChr);
    }

    public static int chromosomeRank(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
            return -1;

        if(chromosome.equals("X"))
            return 23;
        else if(chromosome.equals("Y"))
            return 24;
        else if(chromosome.equals("MT"))
            return 25;
        else
            return Integer.parseInt(chromosome);
    }

}
