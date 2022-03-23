package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_CN_INFO;

import java.util.Set;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class GermlineRescueLowVAF
{
    private static final double MIN_VCN = 0.5;
    private static final int MIN_ALLELE_READ_COUNT = 3;

    private final String mGermlineSample;

    public GermlineRescueLowVAF(final String germlineSample)
    {
        mGermlineSample = germlineSample;
    }

    public VariantContext processVariant(final VariantContext context)
    {
        Set<String> filters = context.getFilters();
        if(filters.size() == 1 && filters.contains(GermlineGenotypeEnrichment.LOW_VAF_FILTER))
        {
            Genotype germlineGenotype = context.getGenotype(mGermlineSample);
            AllelicDepth germlineDepth = AllelicDepth.fromGenotype(germlineGenotype);
            double variantCopyNumber = context.getAttributeAsDouble(PURPLE_VARIANT_CN_INFO, 0.0);
            if(germlineDepth.alleleReadCount() >= MIN_ALLELE_READ_COUNT && Doubles.greaterOrEqual(variantCopyNumber, MIN_VCN))
            {
                return new VariantContextBuilder(context).passFilters().make();
            }
        }

        return context;
    }
}
