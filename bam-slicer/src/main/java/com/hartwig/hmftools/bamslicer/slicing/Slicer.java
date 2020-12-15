package com.hartwig.hmftools.bamslicer.slicing;

import java.util.Collection;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface Slicer extends Predicate<GenomePosition> {

    @NotNull
    Collection<? extends GenomeRegion> regions();

    default boolean includes(@NotNull GenomePosition variant) {
        return test(variant);
    }
}
