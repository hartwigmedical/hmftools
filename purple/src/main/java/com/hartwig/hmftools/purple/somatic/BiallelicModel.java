package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_BASE_ERROR_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_GROWTH_RATE;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_LOH_MIN_PURITY;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_THRESHOLD_PARAMETER_I;
import static com.hartwig.hmftools.purple.PurpleConstants.BIALLELIC_THRESHOLD_PARAMETER_II;

import com.hartwig.hmftools.common.purple.PurpleCopyNumber;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class BiallelicModel
{
    public static double calculateBiallelicProbability(
            final double samplePurity, final PurpleCopyNumber purpleCopyNumber, final SomaticVariant variant)
    {
        int alleleReadCount = variant.alleleReadCount();

        if(alleleReadCount == 0)
            return 0;

        double vaf = variant.alleleFrequency();

        // double minorAlleleCopyNumber = purpleCopyNumber.minorAlleleCopyNumber();
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double variantCopyNumber = variant.decorator().variantCopyNumber();

        double probabilityLoh = probabilityLoh(samplePurity);
        double probabilityNoLoh = 1 - probabilityLoh;

        double vcnThresholdForNoWildtype = calcVcnThresholdForNoWildtype(copyNumber);

        double readCountAtThreshold = calcReadCountAtThreshold(vcnThresholdForNoWildtype, variantCopyNumber, alleleReadCount);

        // dampeningFactor = 1 if alleleReadCount < 300, else (300 + (alleleReadCount - 300)^0.75) / alleleReadCount.

        double probNoWildtypeAssumeLoh = conditionalProbNoWildtypeAssumeLoh(readCountAtThreshold, alleleReadCount);

        double probNoWildtypeAssumeNoLoh = conditionalProbNoWildtypeAssumeNoLoh(probNoWildtypeAssumeLoh, probabilityLoh);

        // Final calculation
        double probabilityNoWildtype = probabilityNoWildtype(
                probabilityLoh, probabilityNoLoh, probNoWildtypeAssumeLoh, probNoWildtypeAssumeNoLoh);

        return probabilityNoWildtype;
    }

    public static boolean classifyBiallelic(double bialleicProbability) { return bialleicProbability >= 0.5; }

    private static double probabilityLoh(double samplePurity)
    {
        return max(samplePurity, BIALLELIC_LOH_MIN_PURITY) * BIALLELIC_LOH_GROWTH_RATE;
    }

    private static double calcVcnThresholdForNoWildtype(double copyNumber)
    {
        double vcnThresholdForNoWildtype = min(copyNumber - 0.5, max(BIALLELIC_THRESHOLD_PARAMETER_I, copyNumber - BIALLELIC_THRESHOLD_PARAMETER_II));
        return vcnThresholdForNoWildtype;
    }

    private static double calcReadCountAtThreshold(double threshold, double variantCopyNumber, int alleleReadCount)
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

    private static double conditionalProbNoWildtypeAssumeNoLoh(double probNoWildtypeAssumeLoh, double probabilityLoh)
    {
        double conditionalProbNoWildtypeAssumeNoLOH = max(probabilityLoh, BIALLELIC_LOH_BASE_ERROR_RATE) / ((1 - probNoWildtypeAssumeLoh)
                + max(probabilityLoh, BIALLELIC_LOH_BASE_ERROR_RATE));

        if(Double.isNaN(conditionalProbNoWildtypeAssumeNoLOH))
        {
            return 0.0d;
        }

        return conditionalProbNoWildtypeAssumeNoLOH;
    }

    private static double probabilityNoWildtype(
            double probabilityLoh, double probabilityNoLoh, double probNoWildtypeAssumeLoh, double probNoWildtypeAssumeNoLoh)
    {
        double probabilityNoWildtype = probabilityLoh * probNoWildtypeAssumeLoh + probabilityNoLoh * probNoWildtypeAssumeNoLoh;

        return probabilityNoWildtype;
    }
}
