package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextEnrichmentComplete implements VariantContextEnrichment {

    public static VariantContextEnrichmentFactory factory(@NotNull final IndexedFastaSequenceFile reference,
            @NotNull final Multimap<String, GenomeRegion> highConfidenceRegions) {
        return consumer -> new VariantContextEnrichmentComplete(reference, highConfidenceRegions, consumer);
    }

    private final KataegisEnrichment kataegisEnrichment;
    private final RefContextEnrichment refContextEnrichment;
    private final HighConfidenceEnrichment highConfidenceEnrichment;

    VariantContextEnrichmentComplete(@NotNull final IndexedFastaSequenceFile reference,
            @NotNull final Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull final Consumer<VariantContext> consumer) {
        highConfidenceEnrichment = new HighConfidenceEnrichment(highConfidenceRegions, consumer);
        kataegisEnrichment = new KataegisEnrichment(highConfidenceEnrichment);
        refContextEnrichment = new RefContextEnrichment(reference, kataegisEnrichment);
    }

    @Override
    public void flush() {
        refContextEnrichment.flush();
        kataegisEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        return kataegisEnrichment.enrichHeader(refContextEnrichment.enrichHeader(template));
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        refContextEnrichment.accept(context);
    }
}
