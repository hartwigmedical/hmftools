package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public interface VariantContextEnrichmentFactory {

    @NotNull
    VariantContextEnrichment create(@NotNull final Consumer<VariantContext> consumer);


    static VariantContextEnrichmentFactory none() {
        return consumer -> new VariantContextEnrichment() {
            @Override
            public void flush() {
                // empty
            }

            @NotNull
            @Override
            public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
                return template;
            }

            @Override
            public void accept(final VariantContext context) {
                consumer.accept(context);
            }
        };
    }

}
