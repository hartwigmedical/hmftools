package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariantHeader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class PurityEnrichment implements VariantContextEnrichment {

    private final String purpleVersion;
    private final PurityAdjustedSomaticVariantFactory factory;
    private final Consumer<VariantContext> consumer;

    PurityEnrichment(@NotNull final String purpleVersion, @NotNull final String tumorSample,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions, @NotNull final Consumer<VariantContext> consumer) {
        this.purpleVersion = purpleVersion;
        this.consumer = consumer;
        this.factory = new PurityAdjustedSomaticVariantFactory(tumorSample, purityAdjuster, copyNumbers, fittedRegions);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        consumer.accept(factory.enrich(context));
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        return SomaticVariantHeader.generateHeader(purpleVersion, template);
    }

    @Override
    public void flush() {
        // Empty
    }
}
