package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_PROB;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_BASE_ERROR_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_GROWTH_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_THRESHOLD_PARAMETER_I;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_THRESHOLD_PARAMETER_II;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.apache.commons.math3.distribution.PoissonDistribution;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticPurityEnrichment
{
    private final PurityAdjuster mPurityAdjuster;
    private final GenomeRegionSelector<PurpleCopyNumber> mCopyNumberSelector;
    private final GenomeRegionSelector<ObservedRegion> mObservedRegionSelector;

    public SomaticPurityEnrichment(
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<ObservedRegion> fittedRegions)
    {
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
        
        double biallelicProbability = calculateBiallelic(purpleCopyNumber, variant);
        boolean classifyBiallelic = classifyBiallelic(biallelicProbability);
        
        PPL_LOGGER.trace("variant({}) biallelic({} prob={})", variant, classifyBiallelic, format("%.4f", biallelicProbability));
        
        variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_PROB, format("%.4f", biallelicProbability)); 
        variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, classifyBiallelic);
    }

    // version 6.0 - New biallelic model
    private static double probabilityLoh(double minorAlleleCopyNumber)
    {
        double probabilityLoh = 1 - 1 / (1 + exp(-BIALLELIC_LOH_GROWTH_RATE * (minorAlleleCopyNumber - 0.5)));

        return probabilityLoh;
    }

    private static double probabilityNoLoh(double probabilityLoh)
    {
        double probabilityNoLoh = 1 - probabilityLoh;

        return probabilityNoLoh;
    }

    private static double vcnThresholdForNoWildtype(double copyNumber)
    {
        double vcnThresholdForNoWildtype = min(copyNumber - 0.5, max(BIALLELIC_THRESHOLD_PARAMETER_I, copyNumber - BIALLELIC_THRESHOLD_PARAMETER_II));

        return vcnThresholdForNoWildtype;
    }

    private static double readCountAtThreshold(double threshold, double variantCopyNumber, int alleleReadCount)
    {
        double readCountAtThreshold = (threshold / variantCopyNumber) * alleleReadCount;

        return readCountAtThreshold;
    }

    private static double conditionalProbNoWildtypeAssumeLoh(double readCountAtThreshold, int alleleReadCount)
    {
        PoissonDistribution poissonDist = new PoissonDistribution(alleleReadCount);

        int readCountAtThresholdInteger = (int) floor(readCountAtThreshold);
        double conditionalProbNoWildtypeAssumeLoh = 1 - poissonDist.cumulativeProbability(readCountAtThresholdInteger);

        return conditionalProbNoWildtypeAssumeLoh;
    }

    private static double conditionalProbNoWildtypeAssumeNoLoh(double conditionalProbNoWildtypeAssumeLoh, double probabilityLoh)
    {
        double conditionalProbNoWildtypeAssumeNoLOH =
                max(probabilityLoh, BIALLELIC_LOH_BASE_ERROR_RATE) / ((1 - conditionalProbNoWildtypeAssumeLoh)
                        + max(probabilityLoh, BIALLELIC_LOH_BASE_ERROR_RATE));

        if(Double.isNaN(conditionalProbNoWildtypeAssumeNoLOH))
        {
            return 0.0d;
        }

        return conditionalProbNoWildtypeAssumeNoLOH;
    }

    private static double probabilityNoWildtype(double probabilityLoh, double probabilityNoLoh, double conditionalProbNoWildtypeAssumeLoh,
            double conditionalProbNoWildtypeAssumeNoLoh)
    {
        double probabilityNoWildtype =
                probabilityLoh * conditionalProbNoWildtypeAssumeLoh + probabilityNoLoh * conditionalProbNoWildtypeAssumeNoLoh;

        return probabilityNoWildtype;
    }

    public static double calculateBiallelic(final PurpleCopyNumber purpleCopyNumber, final SomaticVariant variant)
    {
        // inputs
        double minorAlleleCopyNumber = purpleCopyNumber.minorAlleleCopyNumber();
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double variantCopyNumber = variant.decorator().variantCopyNumber();
        int alleleReadCount = variant.alleleReadCount();

        // part 1
        double probabilityLoh = probabilityLoh(minorAlleleCopyNumber);
        double probabilityNoLoh = probabilityNoLoh(probabilityLoh);

        // part 2
        double vcnThresholdForNoWildtype = vcnThresholdForNoWildtype(copyNumber);
        double readCountAtThreshold = readCountAtThreshold(vcnThresholdForNoWildtype, variantCopyNumber, alleleReadCount);

        // part 3
        double conditionalProbNoWildtypeAssumeLoh = conditionalProbNoWildtypeAssumeLoh(readCountAtThreshold, alleleReadCount);
        double conditionalProbNoWildtypeAssumeNoLoh =
                conditionalProbNoWildtypeAssumeNoLoh(conditionalProbNoWildtypeAssumeLoh, probabilityLoh);

        // Final calculation
        double probabilityNoWildtype =
                probabilityNoWildtype(probabilityLoh, probabilityNoLoh, conditionalProbNoWildtypeAssumeLoh, conditionalProbNoWildtypeAssumeNoLoh);
        return probabilityNoWildtype;
    }

    public static boolean classifyBiallelic(double bialleicProbability) { return bialleicProbability >= 0.5; }
}
