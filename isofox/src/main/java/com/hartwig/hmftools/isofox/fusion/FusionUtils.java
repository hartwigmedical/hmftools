package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class FusionUtils
{

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

    public static String formLocationPair(final String[] chromosomes, final int[] geneCollectionIds)
    {
        return String.format("%s_%s",
                formLocation(chromosomes[SE_START], geneCollectionIds[SE_START]),
                formLocation(chromosomes[SE_END], geneCollectionIds[SE_END]));
    }

    public static String formLocation(final String chromosome, final int geneCollectionId)
    {
        return String.format("%s:%d", chromosome, geneCollectionId);
    }

    public static String formChromosomePair(final String chr1, final String chr2) { return chr1 + "_" + chr2; }
    public static String[] getChromosomePair(final String chrPair) { return chrPair.split("_"); }

}
