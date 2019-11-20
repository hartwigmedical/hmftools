package com.hartwig.hmftools.sage.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public interface SageVariantContextFactory {
    @NotNull
    VariantContext create(@NotNull final SageVariant entry);
}
