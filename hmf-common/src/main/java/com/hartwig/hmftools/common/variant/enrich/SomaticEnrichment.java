package com.hartwig.hmftools.common.variant.enrich;

import java.io.IOException;

import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public interface SomaticEnrichment {

    @NotNull
    default ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull VariantContext context) throws IOException {
        return builder;
    }
}
