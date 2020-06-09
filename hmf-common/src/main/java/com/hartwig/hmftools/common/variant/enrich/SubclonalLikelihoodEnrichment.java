package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariantHeader;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihood;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment implements VariantContextEnrichment {

    public static final String SUBCLONAL_LIKELIHOOD_FLAG = "SUBCL";
    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION = "Non-zero subclonal likelihood";

    private final double maxPloidy;
    private final Consumer<VariantContext> consumer;
    private final SubclonalLikelihood likelihoodFactory;

    SubclonalLikelihoodEnrichment(double maxPloidy, double binWidth, @NotNull final List<PeakModel> peakModel,
            @NotNull final Consumer<VariantContext> consumer) {
        this.maxPloidy = maxPloidy;
        this.consumer = consumer;

        likelihoodFactory = new SubclonalLikelihood(binWidth, peakModel);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        final double variantCopyNumber = context.getAttributeAsDouble(SomaticVariantHeader.PURPLE_VARIANT_CN_INFO, maxPloidy);
        double subclonalLikelihood = Math.round(likelihoodFactory.subclonalLikelihood(variantCopyNumber) * 1000d) / 1000d;

        if (!Doubles.isZero(subclonalLikelihood)) {
            context.getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, subclonalLikelihood);
        }
        consumer.accept(context);
    }

    @Override
    public void flush() {
        // None
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,
                1,
                VCFHeaderLineType.Float, SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION));

        return template;
    }
}
