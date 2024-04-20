package com.hartwig.hmftools.purple.germline;

import com.hartwig.hmftools.common.utils.Doubles;

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

    public static void processVariant(final GermlineVariant variant)
    {
        if(Doubles.lessThan(variant.copyNumber(), MIN_TUMOR_VCN))
        {
            double minQual = variant.isHotspot() ? MIN_QUAL_HOTSPOT : MIN_QUAL_OTHER;

            if(variant.context().getPhredScaledQual() < minQual)
            {
                variant.filters().add(LOW_TUMOR_VCN_FILTER);
                VariantContext newContext = new VariantContextBuilder(variant.context()).filter(LOW_TUMOR_VCN_FILTER).make();
                variant.setContext(newContext);
            }
        }
    }

    public static void enrichHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_TUMOR_VCN_FILTER, "Germline variant has very low tumor variant copy number"));
    }
}
