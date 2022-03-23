package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_CN_INFO;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.Hotspot;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public final class GermlineLowTumorVCNFilter
{
    public static final String LOW_TUMOR_VCN_FILTER = "LOW_TUMOR_VCN";

    public static final double MIN_TUMOR_VCN = 0.5;
    public static final double MIN_QUAL_HOTSPOT = 120;
    public static final double MIN_QUAL_OTHER = 200;

    public static VariantContext processVariant(final VariantContext context)
    {
        double variantCopyNumber = context.getAttributeAsDouble(PURPLE_VARIANT_CN_INFO, 0.0);
        if(Doubles.lessThan(variantCopyNumber, MIN_TUMOR_VCN))
        {
            boolean isHotspot = Hotspot.fromVariant(context) == Hotspot.HOTSPOT;
            double minQual = isHotspot ? MIN_QUAL_HOTSPOT : MIN_QUAL_OTHER;
            if(context.getPhredScaledQual() < minQual)
            {
                return new VariantContextBuilder(context).filter(LOW_TUMOR_VCN_FILTER).make();
            }
        }

        return context;
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFFilterHeaderLine(LOW_TUMOR_VCN_FILTER, "Germline variant has very low tumor variant copy number"));
        return template;
    }
}
