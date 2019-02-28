package com.hartwig.hmftools.common.amber.qc;

import java.util.List;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public class AmberQCFactory {

    @NotNull
    public static AmberQC create(double contamination, @NotNull final List<AmberBAF> baf) {

        final double meanBaf = baf.stream()
                .filter(x -> HumanChromosome.contains(x.chromosome()))
                .filter(x -> HumanChromosome.fromString(x.chromosome()).isAutosome())
                .mapToDouble(AmberBAF::tumorBAF)
                .filter(x -> !Double.isNaN(x))
                .average()
                .orElse(0);

        return ImmutableAmberQC.builder().meanBAF(meanBaf).contamination(contamination).build();
    }
}
