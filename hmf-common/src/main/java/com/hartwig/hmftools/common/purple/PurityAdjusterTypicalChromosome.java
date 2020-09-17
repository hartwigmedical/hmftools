package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;

import org.jetbrains.annotations.NotNull;

public class PurityAdjusterTypicalChromosome extends PurityAdjuster {

    @NotNull
    private final Gender gender;

    public PurityAdjusterTypicalChromosome(@NotNull final Gender gender, @NotNull final FittedPurity fittedPurity) {
        this(gender, fittedPurity.purity(), fittedPurity.normFactor());
    }

    public PurityAdjusterTypicalChromosome(@NotNull final Gender gender, final double purity, final double normFactor) {
        super(purity, normFactor);
        this.gender = gender;
    }

    @Override
    public double germlineRatio(@NotNull final String contig) {
        return HumanChromosome.fromString(contig).isDiploid(gender) ? 1 : 0.5;
    }

}
