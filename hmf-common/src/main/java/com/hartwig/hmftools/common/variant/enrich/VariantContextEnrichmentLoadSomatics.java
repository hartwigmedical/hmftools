package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextEnrichmentLoadSomatics implements VariantContextEnrichment {

    public static VariantContextEnrichmentFactory factory(@NotNull final Multimap<String, GenomeRegion> highConfidenceRegions) {
        return consumer -> new VariantContextEnrichmentLoadSomatics(highConfidenceRegions, consumer);
    }

    private final HighConfidenceEnrichment highConfidenceEnrichment;

    private VariantContextEnrichmentLoadSomatics(@NotNull final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull final Consumer<VariantContext> consumer) {
        highConfidenceEnrichment = new HighConfidenceEnrichment(highConfidenceRegions, consumer);
    }

    @Override
    public void flush() {
        highConfidenceEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        return highConfidenceEnrichment.enrichHeader(template);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        highConfidenceEnrichment.accept(context);
    }
}
