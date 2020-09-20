package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltChromosome extends Chromosome {

    String contig();

    double typicalRatio();

    double actualRatio();

    boolean mosiac();

    @Override
    default boolean isDiploid(@NotNull Gender gender) {
        return Doubles.equal(typicalRatio(), 1.0);
    }

    default boolean isDiploid() {
        return Doubles.equal(typicalRatio(), 1.0);
    }
}
