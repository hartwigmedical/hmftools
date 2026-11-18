package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_PROB;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.Multimaps;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticPurityEnrichment
{
    private final double mSamplePurity;
    private final PurityAdjuster mPurityAdjuster;
    private final GenomeRegionSelector<PurpleCopyNumber> mCopyNumberSelector;
    private final GenomeRegionSelector<ObservedRegion> mObservedRegionSelector;

    public SomaticPurityEnrichment(
            final double samplePurity, final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers,
            final List<ObservedRegion> fittedRegions)
    {
        mSamplePurity = samplePurity;
        mPurityAdjuster = purityAdjuster;
        mCopyNumberSelector = GenomeRegionSelectorFactory.createImproved(Multimaps.fromRegions(copyNumbers));
        mObservedRegionSelector = GenomeRegionSelectorFactory.createImproved(Multimaps.fromRegions(fittedRegions));
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(!HumanChromosome.contains(variant.chromosome()))
            return;

        Optional<ObservedRegion> observedRegion = mObservedRegionSelector.select(variant);
        GermlineStatus germlineStatus = GermlineStatus.UNKNOWN;

        if(observedRegion.isPresent())
            germlineStatus = observedRegion.get().germlineStatus();

        variant.context().getCommonInfo().putAttribute(PURPLE_GERMLINE_INFO, germlineStatus.toString());

        if(variant.hasTumorAlleleDepth())
        {
            Optional<PurpleCopyNumber> purpleCopyNumber = mCopyNumberSelector.select(variant);
            if(purpleCopyNumber.isPresent())
            {
                applyPurityAdjustment(variant, purpleCopyNumber.get(), germlineStatus == GermlineStatus.HET_DELETION);
            }
        }
    }

    private void applyPurityAdjustment(final SomaticVariant variant, final PurpleCopyNumber purpleCopyNumber, boolean isGermlineHetDeletion)
    {
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();

        double vaf = mPurityAdjuster.purityAdjustedVAF(
                purpleCopyNumber.chromosome(), max(0.001, copyNumber), variant.alleleFrequency(), isGermlineHetDeletion);

        double variantCopyNumber = max(0, vaf * copyNumber);
        
        VariantContext variantContext = variant.context();

        variantContext.getCommonInfo().putAttribute(PURPLE_VARIANT_CN, variantCopyNumber);
        variantContext.getCommonInfo().putAttribute(PURPLE_CN, copyNumber);
        variantContext.getCommonInfo().putAttribute(PURPLE_AF, format("%.4f", vaf));
        variantContext.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, purpleCopyNumber.minorAlleleCopyNumber());

        double biallelicProbability = BiallelicModel.calculateProbability(mSamplePurity, purpleCopyNumber, variant);
        boolean classifyBiallelic = BiallelicModel.classifyIsBiallelic(biallelicProbability);

        variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_PROB, format("%.4f", biallelicProbability));
        variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, classifyBiallelic);
    }
}
