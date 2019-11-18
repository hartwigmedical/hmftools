package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public interface VariantContextEnrichmentFactory {

    @NotNull
    VariantContextEnrichment create(@NotNull final Consumer<VariantContext> consumer);

    @NotNull
    static VariantContextEnrichmentFactory noEnrichment() {
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
            public void accept(@NotNull final VariantContext context) {
                consumer.accept(context);
            }
        };
    }
}
