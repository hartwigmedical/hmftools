package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextEnrichmentPurple implements VariantContextEnrichment {

    private final PurityEnrichment purityEnrichment;
    private final KataegisEnrichment kataegisEnrichment;
    private final RefContextEnrichment refContextEnrichment;

    public VariantContextEnrichmentPurple(@NotNull final String purpleVersion, @NotNull final String tumorSample,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<FittedRegion> fittedRegions,
            @NotNull final Consumer<VariantContext> consumer) {
        purityEnrichment = new PurityEnrichment(purpleVersion, tumorSample, purityAdjuster, copyNumbers, fittedRegions, consumer);
        kataegisEnrichment = new KataegisEnrichment(purityEnrichment);
        refContextEnrichment = new RefContextEnrichment(reference, kataegisEnrichment);
    }

    @Override
    public void flush() {
        refContextEnrichment.flush();
        kataegisEnrichment.flush();
        purityEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        VCFHeader header = refContextEnrichment.enrichHeader(template);
        header = kataegisEnrichment.enrichHeader(header);
        return purityEnrichment.enrichHeader(header);
    }

    @Override
    public void accept(final VariantContext context) {
        refContextEnrichment.accept(context);
    }
}
