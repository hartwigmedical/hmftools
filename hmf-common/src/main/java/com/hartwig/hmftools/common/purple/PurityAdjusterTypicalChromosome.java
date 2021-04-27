package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public class PurityAdjusterTypicalChromosome extends PurityAdjuster {

    @NotNull
    private final Gender gender;

    public PurityAdjusterTypicalChromosome(@NotNull final Gender gender, final double purity, final double normFactor) {
        super(purity, normFactor);
        this.gender = gender;
    }

    @Override
    public double germlineRatio(@NotNull final String contig) {
        return HumanChromosome.fromString(contig).isDiploid(gender) ? 1 : 0.5;
    }

}
