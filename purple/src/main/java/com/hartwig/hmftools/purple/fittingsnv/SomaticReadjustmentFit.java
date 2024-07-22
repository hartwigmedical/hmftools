package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.PurpleCopyNumber.buildChromosomeMap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.INVALID_PURITY;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_CN_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_CN_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_PROB_THRESHOLD;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_MINOR_ALLELE_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_READJUST_PURITY_INCREMENT;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class SomaticReadjustmentFit
{
    public static double calcReadjustmentPurity(
            final List<SomaticVariant> variants, final FittedPurity standardSomaticPurity, final List<PurpleCopyNumber> copyNumbers)
    {
        Map<String,List<PurpleCopyNumber>> chrCopyNumberMap = buildChromosomeMap(copyNumbers);

        double standardPurity = standardSomaticPurity.purity();

        List<SomaticVariant> candidateVariants = Lists.newArrayList();
        double maxObservedVaf = 0;

        for(SomaticVariant variant : variants)
        {
            if(!HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            if(!variantInValidCopyNumberRegion(variant, chrCopyNumberMap))
                continue;

            candidateVariants.add(variant);

            maxObservedVaf = max(maxObservedVaf, variant.alleleFrequency());
        }

        if(candidateVariants.size() < SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS)
            return INVALID_PURITY;

        double standardVaf = standardPurity * 0.5;
        List<Double> vafProbabilities = calcVafProbabilities(candidateVariants, standardVaf);

        double probabilityThreshold = calcProbabilityThreshold(vafProbabilities, SNV_READJUST_PROB_THRESHOLD);

        int outlierVafCount = (int)candidateVariants.stream().filter(x -> x.alleleFrequency() >= standardVaf).count();

        double expectedOutlierCount = outlierVafCount * probabilityThreshold / 0.5;

        int outlierCount = (int)vafProbabilities.stream().filter(x -> x >= probabilityThreshold).count();

        int minOutliers = standardPurity < SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL ?
                SNV_READJUST_EXPECTED_VARIANT_COUNT : SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY;

        double requiredCount = max(expectedOutlierCount * SNV_READJUST_EXPECTED_VARIANT_RATIO, minOutliers);

        if(outlierCount < requiredCount)
            return INVALID_PURITY;

        double closestPurity = findClosestExpectedPurity(candidateVariants, standardPurity, maxObservedVaf);

        PPL_LOGGER.debug("reassessment purity({}) from candidate({})", format("%.4f", closestPurity), candidateVariants.size());

        return closestPurity;
    }

    private static boolean variantInValidCopyNumberRegion(final SomaticVariant variant, Map<String,List<PurpleCopyNumber>> chrCopyNumberMap)
    {
        List<PurpleCopyNumber> copyNumbers = chrCopyNumberMap.get(variant.chromosome());

        if(copyNumbers == null)
            return false;

        PurpleCopyNumber copyNumber = copyNumbers.stream()
                .filter(x -> positionWithin(variant.position(), x.start(), x.end())).findFirst().orElse(null);

        if(copyNumber == null)
            return false;

        return copyNumber.averageTumorCopyNumber() >= SNV_READJUST_CN_MIN
            && copyNumber.averageTumorCopyNumber() <= SNV_READJUST_CN_MAX
            && copyNumber.minorAlleleCopyNumber() >= SNV_READJUST_MINOR_ALLELE_MIN;
    }

    private static List<Double> calcVafProbabilities(final List<SomaticVariant> candidateVariants, double testVaf)
    {
        List<Double> vafProbabilities = Lists.newArrayList();

        for(SomaticVariant variant : candidateVariants)
        {
            double expectedAlleleCount = variant.totalReadCount() * testVaf;

            if(variant.alleleReadCount() > expectedAlleleCount)
            {
                BinomialDistribution binomialDist = new BinomialDistribution(variant.totalReadCount(), testVaf);
                double probability = 1 - binomialDist.cumulativeProbability(variant.alleleReadCount() - 1);

                vafProbabilities.add(probability);
            }
            else
            {
                vafProbabilities.add(1.0);
            }
        }

        Collections.sort(vafProbabilities);

        return vafProbabilities;
    }

    private static double calcProbabilityThreshold(final List<Double> vafProbabilities, double probThresholdMin)
    {
        if(vafProbabilities.size() < SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS)
            return probThresholdMin;

        Collections.sort(vafProbabilities);

        double nthVariantProbability = vafProbabilities.get(SNV_READJUST_PROB_THRESHOLD_MIN_VARIANTS - 1);
        return max(nthVariantProbability, probThresholdMin);
    }

    private static double findClosestExpectedPurity(final List<SomaticVariant> candidateVariants, double standardPurity, double maxObservedVaf)
    {
        double maxImpliedPurity = min(maxObservedVaf * 2, 1);

        // calculate expected number of ‘outliers’ calculated as 1% of the observed variants with VAF > 0.5 * candidate purity
        // Note that if there are <500 variants with VAF > 0.5, the outlier percentile is recalculated as the percentile that would give 5 expected outliers.
        // The sample purity is then readjusted to be the lowest purity which satisfies observed outliers < expected outliers.

        double minValidPurity = INVALID_PURITY;

        for(double testPurity = standardPurity; testPurity < maxImpliedPurity; testPurity += SNV_READJUST_PURITY_INCREMENT)
        {
            double testVaf = testPurity * 0.5;

            List<Double> vafProbabilities = calcVafProbabilities(candidateVariants, testVaf);

            double probabilityThreshold = calcProbabilityThreshold(vafProbabilities, SNV_READJUST_PROB_THRESHOLD);

            int outlierVafCount = (int)candidateVariants.stream().filter(x -> x.alleleFrequency() >= testVaf).count();

            double expectedOutlierCount = outlierVafCount * probabilityThreshold / 0.5;

            int outlierCount = (int)vafProbabilities.stream().filter(x -> x <= probabilityThreshold).count();

            if(outlierCount <= expectedOutlierCount)
            {
                if(minValidPurity == INVALID_PURITY || testPurity < minValidPurity)
                    minValidPurity = testPurity;
            }
        }

        return minValidPurity;
    }
}
