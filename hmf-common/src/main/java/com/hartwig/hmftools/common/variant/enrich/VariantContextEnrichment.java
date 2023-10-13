package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public interface VariantContextEnrichment extends Consumer<VariantContext>
{
    @Override
    void accept(@NotNull final VariantContext context);

    void flush();

    @NotNull
    VCFHeader enrichHeader(@NotNull final VCFHeader template);
}
