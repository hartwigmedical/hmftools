package com.hartwig.hmftools.common.slicing;

import java.util.Collection;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public interface Slicer extends Predicate<GenomePosition> {

    @NotNull
    Collection<? extends GenomeRegion> regions();

    default int numberOfRegions() {
        return regions().size();
    }

    default long numberOfBases() {
        return regions().stream().mapToLong(GenomeRegion::bases).sum();
    }

    default boolean includes(@NotNull Variant variant) {
        return test(variant);
    }
}
