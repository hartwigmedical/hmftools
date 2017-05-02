package com.hartwig.hmftools.common.slicing;

import com.hartwig.hmftools.common.variant.Variant;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.function.Predicate;

public interface Slicer extends Predicate<Variant> {

    @NotNull
    Collection<GenomeRegion> regions();

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
