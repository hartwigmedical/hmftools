package com.hartwig.hmftools.common.variant.enrich;

import java.io.IOException;
import java.util.ArrayList;

import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class CompoundEnrichment extends ArrayList<SomaticEnrichment> implements SomaticEnrichment {

    public CompoundEnrichment() {
        super();
    }

    public CompoundEnrichment(@NotNull SomaticEnrichment initialEnrichment) {
        super();
        add(initialEnrichment);
    }

    @NotNull
    @Override
    public ImmutableSomaticVariantImpl.Builder enrich(@NotNull final ImmutableSomaticVariantImpl.Builder builder,
            @NotNull final VariantContext context) throws IOException {
        for (SomaticEnrichment enrichment : this) {
            enrichment.enrich(builder, context);
        }

        return builder;
    }
}
