package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_AF_UNCERTAINTY_FACTOR;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_DAMPENING_AD_EXP;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_DAMPENING_AD_LIMIT;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_BASE_ERROR_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_GROWTH_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_MIN_PURITY;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_VCN_THRESHOLD_PI;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_VCN_THRESHOLD_P2;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class BiallelicModel
{
    public static boolean classifyBiallelic(double bialleicProbability) { return bialleicProbability >= 0.5; }

    public static double calcBiallelicProbability(
            final double samplePurity, final PurpleCopyNumber purpleCopyNumber, final SomaticVariant variant)
    {
        int alleleReadCount = variant.alleleReadCount();

        if(alleleReadCount == 0)
            return 0;

        int depth = variant.totalReadCount();
        double vaf = variant.alleleFrequency();

        // double minorAlleleCopyNumber = purpleCopyNumber.minorAlleleCopyNumber();
        double segmentCN = purpleCopyNumber.averageTumorCopyNumber();
        double variantCN = variant.decorator().variantCopyNumber();
        double majorAlleleCN = variantCN - variant.decorator().minorAlleleCopyNumber();

        return calcBiallelicProbability(samplePurity, alleleReadCount, depth, vaf, segmentCN, variantCN, majorAlleleCN);
    }

    public static double calcBiallelicProbability(
            final double samplePurity, int alleleReadCount, int depth, double vaf, double segmentCN, double variantCN, double majorAlleleCN)
    {
        double lohExponent = max(samplePurity, BIALLELIC_LOH_MIN_PURITY) * BIALLELIC_LOH_GROWTH_RATE;
        double probabilityLoh = 1 - 1 / (1 + Math.exp(-lohExponent * (majorAlleleCN - 0.5)));
        // =1-1/(1+EXP(-E3*(B19-0.5)))

        double vcnThresholdForNoWildtype = min(segmentCN - 0.5, max(BIALLELIC_VCN_THRESHOLD_PI, segmentCN - BIALLELIC_VCN_THRESHOLD_P2));

        double readCountAtThreshold = vcnThresholdForNoWildtype / variantCN * alleleReadCount;

        double afUncertainty = vaf * BIALLELIC_AF_UNCERTAINTY_FACTOR;

        double conservativeAd = round(Doubles.median(List.of(
                readCountAtThreshold,
                alleleReadCount - depth * afUncertainty,
                alleleReadCount + depth * afUncertainty)));

        if(conservativeAd > BIALLELIC_DAMPENING_AD_LIMIT)
        {
            double factor = pow(conservativeAd - BIALLELIC_DAMPENING_AD_LIMIT, BIALLELIC_DAMPENING_AD_EXP);
            double dampeningFactor = (BIALLELIC_DAMPENING_AD_LIMIT + factor) / conservativeAd;

            conservativeAd = round(conservativeAd * dampeningFactor);
            readCountAtThreshold = round(readCountAtThreshold * dampeningFactor);
        }

        PoissonDistribution poissonDist = new PoissonDistribution(conservativeAd);

        double probNoWildtypeOrLoh = 1 - poissonDist.cumulativeProbability((int)readCountAtThreshold);

        double probWildtypeOrLohRelative = probNoWildtypeOrLoh == 0 ? PROB_UPPER_LIMIT : 0.5 * (1 - probNoWildtypeOrLoh) / probNoWildtypeOrLoh;

        double probabilityLohFloored = max(probabilityLoh, BIALLELIC_LOH_BASE_ERROR_RATE);
        double probNoWildtypeOrNoLohRaw = probabilityLohFloored / (probWildtypeOrLohRelative + probabilityLohFloored);

        double minAlleleDiff = max(min(abs(variantCN - majorAlleleCN), abs(variantCN - (segmentCN - majorAlleleCN))), 1 / PROB_UPPER_LIMIT);
        double minCNDiff = max(abs(variantCN - segmentCN), 1 / PROB_UPPER_LIMIT);

        double probNoWildtypeOrNoLohAdjusted = Math.log(minAlleleDiff / minCNDiff);

        double probNoWildtypeOrNoLoh = calcScaledAdjustment(probNoWildtypeOrNoLohRaw, probNoWildtypeOrNoLohAdjusted);
        probNoWildtypeOrNoLoh = min(max(probNoWildtypeOrNoLoh, 0), 1);

        double biallelicProb = probabilityLoh * probNoWildtypeOrLoh + (1 - probabilityLoh) * probNoWildtypeOrNoLoh;

        return biallelicProb;
    }

    private static final double PROB_UPPER_LIMIT = 1e10;

    private static double calcScaledAdjustment(double origValue, double adjustment, boolean upper)
    {
        return upper ? 1 - (1 - origValue) / pow(Math.E, adjustment) : origValue * pow(Math.E, adjustment);
    }

    private static double calcScaledAdjustment(double origValue, double origAdjustment)
    {
        // propensity to move is inversely proportional to distance to (0, 1) as values near 0 and 1 are already confident
        // conceptually, imagine 1/x for x < 0.5 and 1/(1-x) for x > 0.5. need to find the point where the integral of the change in
        // orig reaches the adjustment

        if(origValue > 0.5)
        {
            if(origAdjustment > 0)
            {
                return calcScaledAdjustment(origValue, origAdjustment, true);
            }
            else
            {
                double areaToMidpoint = -log(0.5) + log(1 - origValue);

                if(origAdjustment - areaToMidpoint > 0)
                {
                    return calcScaledAdjustment(origValue, origAdjustment, true);
                }
                else
                {
                    return calcScaledAdjustment(0.5, origAdjustment - areaToMidpoint, false);
                }
            }
        }
        else
        {
            if(origAdjustment < 0)
            {
                return calcScaledAdjustment(origValue, origAdjustment, false);
            }
            else
            {
                double areaToMidpoint = log(0.5) - log(origValue);

                if(origAdjustment - areaToMidpoint < 0)
                {
                    return calcScaledAdjustment(origValue, origAdjustment, false);
                }
                else
                {
                    return calcScaledAdjustment(0.5, origAdjustment - areaToMidpoint, true);
                }
            }
        }
    }
}
