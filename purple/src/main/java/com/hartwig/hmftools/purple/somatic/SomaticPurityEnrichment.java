package com.hartwig.hmftools.purple.somatic;

import java.util.List;

import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.variant.VariantHeader;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticPurityEnrichment
{
    private final String mPurpleVersion;
    private final PurityAdjustedSomaticVariantFactory mFactory;

    public SomaticPurityEnrichment(final String purpleVersion, final String sample, final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
        mPurpleVersion = purpleVersion;
        mFactory = new PurityAdjustedSomaticVariantFactory(sample, purityAdjuster, copyNumbers, fittedRegions);
    }

    public void processVariant(final VariantContext context)
    {
        mFactory.enrich(context);
    }

    public VCFHeader enrichHeader(final VCFHeader template) {
        return VariantHeader.somaticHeader(mPurpleVersion, template);
    }
}
