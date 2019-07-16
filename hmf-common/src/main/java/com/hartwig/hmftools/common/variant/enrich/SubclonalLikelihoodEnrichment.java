package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihood;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihoodFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment implements VariantContextEnrichment {

    public static final String SUBCLONAL_LIKELIHOOD_FLAG = "SUBCL";
    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCIPTION = "Subclonal likelihood";

    private final double maxPloidy;
    private final double binWidth;
    private final Consumer<VariantContext> consumer;
    private final double[] subclonalLikelihoods;

    SubclonalLikelihoodEnrichment(double maxPloidy, double binWidth, @NotNull final List<PeakModel> peakModel,
            @NotNull final Consumer<VariantContext> consumer) {
        this.maxPloidy = maxPloidy;
        this.binWidth = binWidth;
        this.consumer = consumer;
        this.subclonalLikelihoods = new double[(int) Math.round(maxPloidy / binWidth)];

        final List<SubclonalLikelihood> likelihoods = SubclonalLikelihoodFactory.subclonalLikelihood(binWidth, peakModel);
        for (SubclonalLikelihood likelihood : likelihoods) {
            if (!Doubles.isZero(likelihood.likelihood())) {
                int key = (int) Math.round(likelihood.bucket() / binWidth);
                subclonalLikelihoods[key] = likelihood.likelihood();
            }
        }
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        final double ploidy = context.getAttributeAsDouble(PurityEnrichment.PURPLE_PLOIDY_INFO, maxPloidy);
        int key = (int) Math.round(ploidy / binWidth);
        if (key < subclonalLikelihoods.length) {
            double variantLikelihood = Math.round(subclonalLikelihoods[key] * 1000d) / 1000d;
            if (!Doubles.isZero(variantLikelihood)) {
                context.getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, variantLikelihood);
            }
        }
        consumer.accept(context);
    }

    @Override
    public void flush() {

    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,
                1,
                VCFHeaderLineType.Float,
                SUBCLONAL_LIKELIHOOD_FLAG_DESCIPTION));

        return template;
    }
}
