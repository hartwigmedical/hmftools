package com.hartwig.hmftools.purple.sv;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public interface VariantContextCollection extends Iterable<VariantContext> {

    void add(@NotNull final VariantContext variantContext);

    int remove(@NotNull final Predicate<VariantContext> removePredicate);

    @NotNull
    List<StructuralVariant> segmentationVariants();
}
