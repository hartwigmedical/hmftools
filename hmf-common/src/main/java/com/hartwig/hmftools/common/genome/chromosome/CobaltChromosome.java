package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltChromosome extends Chromosome
{
    @NotNull
    String contig();

    double typicalRatio();

    double actualRatio();

    boolean mosiac();

    default boolean isNormal() {
        return Doubles.equal(typicalRatio(), actualRatio());
    }

    default boolean isDiploid() {
        return Doubles.equal(actualRatio(), 1.0);
    }
}
