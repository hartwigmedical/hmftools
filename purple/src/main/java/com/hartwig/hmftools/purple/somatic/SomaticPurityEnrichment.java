package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantHeader;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticPurityEnrichment
{
    private final String mPurpleVersion;
    private final PurityAdjustedSomaticVariantFactory mFactory;

    public SomaticPurityEnrichment(final String purpleVersion, final String sample, final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions)
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
