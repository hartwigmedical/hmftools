package com.hartwig.hmftools.purple.purity;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;

import org.jetbrains.annotations.NotNull;

public class PurityAdjusterTypicalChromosome extends PurityAdjuster
{
    @NotNull
    private final Gender mGender;

    public PurityAdjusterTypicalChromosome(@NotNull final Gender gender, final double purity, final double normFactor)
    {
        super(purity, normFactor);
        mGender = gender;
    }

    @Override
    public double germlineRatio(@NotNull final String contig)
    {
        return HumanChromosome.fromString(contig).isDiploid(mGender) ? 1 : 0.5;
    }

}
