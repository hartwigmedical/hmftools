package com.hartwig.hmftools.common.variant.enrich;

import java.util.ArrayList;

import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class CompoundEnrichment extends ArrayList<SomaticEnrichment> implements SomaticEnrichment {
    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) {
        for (SomaticEnrichment enrichment : this) {
            enrichment.enrich(builder, context);
        }

        return builder;
    }
}
