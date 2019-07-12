package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextEnrichmentPurple implements VariantContextEnrichment {

    private final PurityEnrichment purityEnrichment;
    private final KataegisEnrichment kataegisEnrichment;
    private final SomaticRefContextEnrichment somaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment subclonalLikelihoodEnrichment;

    public VariantContextEnrichmentPurple(@NotNull final String purpleVersion, @NotNull final String tumorSample,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<FittedRegion> fittedRegions,
            @NotNull final Consumer<VariantContext> consumer) {
        subclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(tumorSample, consumer);
        purityEnrichment =
                new PurityEnrichment(purpleVersion, tumorSample, purityAdjuster, copyNumbers, fittedRegions, subclonalLikelihoodEnrichment);
        kataegisEnrichment = new KataegisEnrichment(purityEnrichment);
        somaticRefContextEnrichment = new SomaticRefContextEnrichment(reference, kataegisEnrichment);
    }

    @NotNull
    public List<PeakModel> clonalityModel() {
        return subclonalLikelihoodEnrichment.peakModel();
    }

    @Override
    public void flush() {
        somaticRefContextEnrichment.flush();
        kataegisEnrichment.flush();
        purityEnrichment.flush();
        subclonalLikelihoodEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        VCFHeader header = somaticRefContextEnrichment.enrichHeader(template);
        header = kataegisEnrichment.enrichHeader(header);
        header = subclonalLikelihoodEnrichment.enrichHeader(header);
        return purityEnrichment.enrichHeader(header);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        somaticRefContextEnrichment.accept(context);
    }
}
