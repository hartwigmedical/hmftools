package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN_INFO;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.PurpleVcfTags;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticPurityEnrichment
{
    private final PurityAdjuster mPurityAdjuster;
    private final GenomeRegionSelector<PurpleCopyNumber> mCopyNumberSelector;
    private final GenomeRegionSelector<ObservedRegion> mObservedRegionSelector;
    private final String mSample;
    private final String mPurpleVersion;

    public SomaticPurityEnrichment(final String purpleVersion, final String sample, final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
        mPurpleVersion = purpleVersion;

        mSample = sample;
        mPurityAdjuster = purityAdjuster;
        mCopyNumberSelector = GenomeRegionSelectorFactory.createImproved(Multimaps.fromRegions(copyNumbers));
        mObservedRegionSelector = GenomeRegionSelectorFactory.createImproved(Multimaps.fromRegions(fittedRegions));
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(variant.hasTumorAlleleDepth())
        {
            Optional<PurpleCopyNumber> purpleCopyNumber = mCopyNumberSelector.select(variant);
            if(purpleCopyNumber.isPresent())
            {
                applyPurityAdjustment(variant, purpleCopyNumber.get());
            }
        }

        Optional<ObservedRegion> observedRegion = mObservedRegionSelector.select(variant);

        if(observedRegion.isPresent())
        {
            variant.context().getCommonInfo().putAttribute(PURPLE_GERMLINE_INFO, observedRegion.get().germlineStatus().toString());
        }
    }

    private void applyPurityAdjustment(final SomaticVariant variant, final PurpleCopyNumber purpleCopyNumber)
    {
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double vaf = mPurityAdjuster.purityAdjustedVAF(purpleCopyNumber.chromosome(), Math.max(0.001, copyNumber), variant.alleleFrequency());
        double ploidy = Math.max(0, vaf * copyNumber);

        boolean biallelic = Doubles.lessOrEqual(copyNumber, 0) || Doubles.greaterOrEqual(ploidy, copyNumber - 0.5);

        VariantContext variantContext = variant.context();

        variantContext.getCommonInfo().putAttribute(PURPLE_VARIANT_CN_INFO, ploidy);
        variantContext.getCommonInfo().putAttribute(PURPLE_CN_INFO, copyNumber);

        variantContext.getCommonInfo().putAttribute(PURPLE_AF_INFO, String.format("%.4f", vaf));
        variantContext.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, purpleCopyNumber.minorAlleleCopyNumber());
        variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, biallelic);
    }

    public VCFHeader enrichHeader(final VCFHeader template)
    {
        return PurpleVcfTags.addSomaticHeader(mPurpleVersion, template);
    }
}
