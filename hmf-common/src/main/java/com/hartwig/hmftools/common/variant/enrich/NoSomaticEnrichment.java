package com.hartwig.hmftools.common.variant.enrich;

import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class NoSomaticEnrichment implements SomaticEnrichment {

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) {
        return builder;
    }
}
