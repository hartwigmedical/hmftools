package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.PurpleCopyNumber.buildChromosomeMap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.INVALID_PURITY;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_CN_MAX;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_CN_MIN;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_EXPECTED_VARIANT_RATIO;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_MINOR_ALLELE_MIN_MIN;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_PROB_THRESHOLD;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_READJUST_PURITY_INCREMENT;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class SomaticReadjustmentFit
{
    public SomaticReadjustmentFit()
    {

    }

    public static double calcReadjustmentPurity(
            final List<SomaticVariant> variants, final FittedPurity standardSomaticPurity, final List<PurpleCopyNumber> copyNumbers)
    {
        Map<String,List<PurpleCopyNumber>> chrCopyNumberMap = buildChromosomeMap(copyNumbers);

        double standardPurity = standardSomaticPurity.purity();

        List<SomaticVariant> candidateVariants = Lists.newArrayList();
        // List<CandidateOutlier> candidateOutliers = Lists.newArrayList();
        // List<Double> expectedAlleleCounts = Lists.newArrayList();
        double maxObservedVaf = 0;

        for(SomaticVariant variant : variants)
        {
            if(!HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            if(!variantInValidCopyNumberRegion(variant, chrCopyNumberMap))
                continue;

            candidateVariants.add(variant);

            maxObservedVaf = max(maxObservedVaf, variant.alleleFrequency());
            /*
            double expectedAlleleCount = variant.totalReadCount() * standardPurity * 0.5;
            expectedAlleleCounts.add(expectedAlleleCount);

            if(variant.alleleReadCount() > expectedAlleleCount)
            {
                PoissonDistribution poissonDist = new PoissonDistribution(expectedAlleleCount);
                double probability = 1 - poissonDist.cumulativeProbability(variant.alleleReadCount() - 1);

                CandidateOutlier candidateOutlier = new CandidateOutlier(variant, probability);
                candidateOutliers.add(candidateOutlier);
            }

            PPL_LOGGER.trace(format("hotspot(%s:%d) vaf(%.3f %d/%d)",
                    variant.chromosome(), variant.position(),
                    variant.alleleFrequency(), variant.alleleReadCount(), variant.totalReadCount()));
            */
        }

        int outlierCount = calcOutlierCount(candidateVariants, standardPurity);

        double expectedOutlierCount = candidateVariants.size() * SNV_READJUST_PROB_THRESHOLD;
        int minOutliers = standardPurity >= SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY_LEVEL ?
                SNV_READJUST_EXPECTED_VARIANT_COUNT : SNV_READJUST_EXPECTED_VARIANT_COUNT_LOW_PURITY;

        double requiredCount = max(expectedOutlierCount * SNV_READJUST_EXPECTED_VARIANT_RATIO, minOutliers);

        /*
        Collections.sort(candidateOutliers, Comparator.comparingDouble(x -> x.Probability));

        // find the nth percentile
        CandidateOutlier nthCandidate = candidateOutliers.get(nthPercIndex);

        PPL_LOGGER.debug("nth outlier candidate({})", nthCandidate.Variant);

        // find the nth percentile of expected counts
        int nthPercIndex = (int)floor(expectedAlleleCounts.size() * SNV_READJUST_NTH_PERC);
        Collections.reverse(expectedAlleleCounts);
        double nthExpectedCount = expectedAlleleCounts.get(nthPercIndex);

        int outlierCount = (int)candidateOutliers.stream().filter(x -> x.Variant.alleleReadCount() > nthExpectedCount).count();
        */

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
            && copyNumber.minorAlleleCopyNumber() >= SNV_READJUST_MINOR_ALLELE_MIN_MIN;
    }

    private static int calcOutlierCount(final List<SomaticVariant> candidateVariants, double purity)
    {
        int outlierCount = 0;

        for(SomaticVariant variant : candidateVariants)
        {
            double expectedAlleleCount = variant.totalReadCount() * purity * 0.5;

            if(variant.alleleReadCount() > expectedAlleleCount)
            {
                PoissonDistribution poissonDist = new PoissonDistribution(expectedAlleleCount);
                double probability = 1 - poissonDist.cumulativeProbability(variant.alleleReadCount() - 1);

                if(probability <= SNV_READJUST_PROB_THRESHOLD)
                    ++outlierCount;
            }
        }

        return outlierCount;
    }

    private static double calcExpectedOutlierCount(final List<SomaticVariant> candidateVariants, double purity)
    {
        double variantAboveVafLevel = candidateVariants.stream().filter(x -> x.alleleFrequency() > purity * 0.5).count();
        return variantAboveVafLevel * 0.01;
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
            int outlierCount = calcOutlierCount(candidateVariants, testPurity);

            double expectedOutlierCount = calcExpectedOutlierCount(candidateVariants, testPurity);

            if(outlierCount <= expectedOutlierCount)
            {
                if(minValidPurity == INVALID_PURITY || testPurity < minValidPurity)
                    minValidPurity = testPurity;
            }
        }

        return minValidPurity;
    }

    private class CandidateOutlier
    {
        public final SomaticVariant Variant;
        public final double Probability;

        public CandidateOutlier(final SomaticVariant variant, final double probability)
        {
            Variant = variant;
            Probability = probability;
        }

        public double vaf() { return Variant.alleleReadCount() / (double)Variant.totalReadCount(); }

        public String toString() { return format("%s prob(%.6f)", Variant, Probability); }
    }
}
