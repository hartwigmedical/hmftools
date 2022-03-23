package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.SUBCLONAL_LIKELIHOOD_FLAG;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.purple.fitting.PeakModel;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment
{
    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION = "Non-zero subclonal likelihood";

    private final SubclonalLikelihood mSubclonalLikelihood;

    public SubclonalLikelihoodEnrichment(double binWidth, final List<PeakModel> peakModel)
    {
        mSubclonalLikelihood = new SubclonalLikelihood(binWidth, peakModel);
    }

    public void processVariant(final VariantContext context)
    {
        VariantContextDecorator decorator = new VariantContextDecorator(context);
            double subclonalLikelihood = Math.round(mSubclonalLikelihood.subclonalLikelihood(decorator.variantCopyNumber()) * 1000d) / 1000d;

        if(!Doubles.isZero(subclonalLikelihood))
        {
            context.getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, subclonalLikelihood);
        }
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(
                new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,1, VCFHeaderLineType.Float, SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION));

        return template;
    }
}
