package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION;

import java.util.List;

import com.hartwig.hmftools.purple.fittingsnv.PeakModelData;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment
{
    private final SubclonalLikelihood mSubclonalLikelihood;

    public SubclonalLikelihoodEnrichment(double binWidth, final List<PeakModelData> peakModel)
    {
        mSubclonalLikelihood = new SubclonalLikelihood(binWidth, peakModel);
    }

    public void processVariant(final SomaticVariant variant)
    {
        double copyNumber = variant.copyNumber();
        double subclonalLikelihood = Math.round(mSubclonalLikelihood.subclonalLikelihood(copyNumber) * 1000d) / 1000d;
        variant.context().getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, subclonalLikelihood);
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(
                new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,1, VCFHeaderLineType.Float, SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION));
    }
}
