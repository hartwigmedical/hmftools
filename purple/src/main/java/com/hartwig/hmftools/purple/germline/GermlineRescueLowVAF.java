package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import htsjdk.variant.variantcontext.Genotype;

public class GermlineRescueLowVAF
{
    private static final double MIN_VCN = 0.5;
    private static final int MIN_ALLELE_READ_COUNT = 3;

    private final String mGermlineSample;

    public GermlineRescueLowVAF(final String germlineSample)
    {
        mGermlineSample = germlineSample;
    }

    public void processVariant(final GermlineVariant variant)
    {
        if(variant.filters().size() == 1 && variant.filters().contains(GermlineGenotypeEnrichment.LOW_VAF_FILTER))
        {
            Genotype germlineGenotype = variant.getGenotype(mGermlineSample);
            AllelicDepth germlineDepth = AllelicDepth.fromGenotype(germlineGenotype);

            if(germlineDepth.AlleleReadCount >= MIN_ALLELE_READ_COUNT && Doubles.greaterOrEqual(variant.copyNumber(), MIN_VCN))
            {
                variant.filters().clear();
                variant.filters().add(PASS);
            }
        }
    }
}
