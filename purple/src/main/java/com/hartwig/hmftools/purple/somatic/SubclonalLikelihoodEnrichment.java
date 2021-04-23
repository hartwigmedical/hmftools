package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.enrich.Subclonality.SUBCLONAL_LIKELIHOOD_FLAG;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment implements VariantContextEnrichment
{

    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION = "Non-zero subclonal likelihood";

    private final Consumer<VariantContext> consumer;
    private final SubclonalLikelihood likelihoodFactory;

    SubclonalLikelihoodEnrichment(double binWidth, final List<PeakModel> peakModel, @NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
        likelihoodFactory = new SubclonalLikelihood(binWidth, peakModel);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        VariantContextDecorator decorator = new VariantContextDecorator(context);
        double subclonalLikelihood = Math.round(likelihoodFactory.subclonalLikelihood(decorator.variantCopyNumber()) * 1000d) / 1000d;

        if (!Doubles.isZero(subclonalLikelihood)) {
            context.getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, subclonalLikelihood);
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
                SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION));

        return template;
    }
}
