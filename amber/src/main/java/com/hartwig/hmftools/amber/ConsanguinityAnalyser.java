package com.hartwig.hmftools.amber;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.Nullable;

public class ConsanguinityAnalyser
{
    // ConsanguinityProportion as the sum of homozygous regions > 3MB divided by the size of the
    // autosome genome (2.88bn bases).
    public static double calcConsanguinityProportion(List<RegionOfHomozygosity> rohs)
    {
        double longRohLengthSum = 0;
        for (RegionOfHomozygosity roh : rohs)
        {
            if (roh.getLength() > AmberConstants.HOMOZYGOUS_REGION_LONG_SIZE)
            {
                longRohLengthSum += roh.getLength();
            }
        }

        return longRohLengthSum / 2_880_000_000L;
    }

    // UniparentalDisomy is if one chromosome has long stretch of homozygosity that
    // means it likely comes from one parent
    // All long roh has to come from one chromosome, and together add up to 10M bases
    @Nullable
    public static Chromosome findUniparentalDisomy(List<RegionOfHomozygosity> rohs)
    {
        Chromosome longRohChromosome = null;
        long longRohLengthSum = 0;

        for (RegionOfHomozygosity roh : rohs)
        {
            if (roh.getLength() >= AmberConstants.HOMOZYGOUS_REGION_LONG_SIZE)
            {
                if (longRohChromosome == null)
                {
                    longRohChromosome = roh.getChromosome();
                }
                else if (!roh.getChromosome().equals(longRohChromosome))
                {
                    // more than 1 chromosome has long ROHs
                    return null;
                }

                longRohLengthSum += roh.getLength();
            }
        }

        if (longRohLengthSum >= AmberConstants.UNIPARENTAL_DISOMY_MIN_LENGTH)
        {
            // if sum of ROH length is long enough
            return longRohChromosome;
        }
        return null;
    }
}
