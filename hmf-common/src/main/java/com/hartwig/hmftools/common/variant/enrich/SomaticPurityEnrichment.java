package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantHeader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticPurityEnrichment implements VariantContextEnrichment {

    private final String purpleVersion;
    private final PurityAdjustedSomaticVariantFactory factory;
    private final Consumer<VariantContext> consumer;

    SomaticPurityEnrichment(@NotNull final String purpleVersion, @NotNull final String sample,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions, @NotNull final Consumer<VariantContext> consumer) {
        this.purpleVersion = purpleVersion;
        this.consumer = consumer;
        this.factory = new PurityAdjustedSomaticVariantFactory(sample, purityAdjuster, copyNumbers, fittedRegions);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        consumer.accept(factory.enrich(context));
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        return VariantHeader.somaticHeader(purpleVersion, template);
    }

    @Override
    public void flush() {
        // Empty
    }
}
