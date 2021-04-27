package com.hartwig.hmftools.purple.fitting;

import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;

import org.jetbrains.annotations.NotNull;

public final class WholeGenomeDuplication
{
    public static final double MIN_AVERAGE_PLOIDY = 1.5;
    public static final int MIN_DUPLICATED_AUTOSOMES = 11;

    public static boolean wholeGenomeDuplication(@NotNull final List<PurpleCopyNumber> copyNumbers)
    {
        return duplicatedAutosomes(copyNumbers) >= MIN_DUPLICATED_AUTOSOMES;
    }

    private static int duplicatedAutosomes(@NotNull final List<PurpleCopyNumber> copyNumbers)
    {
        ListMultimap<Chromosome, PurpleCopyNumber> copyNumberMap = Multimaps.fromRegions(copyNumbers);

        int duplicatedAutosomes = 0;
        for(Chromosome chromosome : copyNumberMap.keySet())
        {
            if(chromosome.isAutosome() && Doubles.greaterOrEqual(averageMajorAlleleCopyNumber(copyNumberMap.get(chromosome)),
                    MIN_AVERAGE_PLOIDY))
            {
                duplicatedAutosomes++;
            }
        }

        return duplicatedAutosomes;
    }

    public static double averageMajorAlleleCopyNumber(@NotNull final List<PurpleCopyNumber> copyNumbers)
    {
        double weightedMajorAllelePloidy = 0;
        long totalBafCount = 0;

        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            weightedMajorAllelePloidy += copyNumber.majorAlleleCopyNumber() * copyNumber.bafCount();
            totalBafCount += copyNumber.bafCount();
        }

        return totalBafCount > 0 ? weightedMajorAllelePloidy / totalBafCount : 0;
    }
}
