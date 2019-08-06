package com.hartwig.hmftools.purple.sv;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextCollectionDummy implements VariantContextCollection {

    @Override
    public void add(@NotNull final VariantContext variantContext) {
    }

    @Override
    public int remove(@NotNull final Predicate<VariantContext> removePredicate) {
        return 0;
    }

    @NotNull
    @Override
    public List<StructuralVariant> segmentationVariants() {
        return Collections.emptyList();
    }

    @NotNull
    @Override
    public Iterator<VariantContext> iterator() {
        return Collections.emptyIterator();
    }
}
