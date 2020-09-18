package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltChromosome extends Chromosome {

    String contig();

    int typicalCopies();

    int impliedCopies();

    boolean mosiac();

    @Override
    default boolean isDiploid(@NotNull Gender gender) {
        return typicalCopies() == 2;
    }
}
