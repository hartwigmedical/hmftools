package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.abs;
import static java.lang.Math.round;
import static java.lang.Math.sin;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.PurpleCommon.formatDbl;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.ImmutableSomaticPeak;
import com.hartwig.hmftools.common.purple.purity.SomaticPeak;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.WeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.WeightedPloidyHistogram;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class SomaticHistogramPeaks
{
    private static final int MAX_ITERATIONS = 10;
    private static final double MAX_UNEXPLAINED_WEIGHT_PERCENT = 0.01;

    private final double mMaxVaf;
    private final double mBinWidth;
    private final double mMinPeakWeight;
    private final double mMinPeakPerc;
    private final int mPeakWidth;

    private final WeightedPloidyHistogram mHistogram;
    private final Map<String, BinomialDistribution> mBinomialDistributionMap;

    public SomaticHistogramPeaks(final double maxVaf, double binWidth, int peakWidth, double minPeakWeight, double minPeakPerc)
    {
        mBinWidth = binWidth;
        mMaxVaf = maxVaf;
        mPeakWidth = peakWidth;
        mMinPeakWeight = minPeakWeight;
        mMinPeakPerc = minPeakPerc;

        mHistogram = new WeightedPloidyHistogram(mMaxVaf, mBinWidth);
        mBinomialDistributionMap = Maps.newHashMap();
    }

    public double findVafPeak(final List<WeightedPloidy> weightedVAFs)
    {
        // find the highest VAF peak
        final List<ModifiableWeightedPloidy> weightedVariants = weightedVAFs.stream()
                .map(x -> ModifiableWeightedPloidy.create()
                        .setTotalReadCount(x.totalReadCount())
                        .setAlleleReadCount(x.alleleReadCount())
                        .setPloidy(x.ploidy())
                        .setWeight(x.weight()))
                .collect(Collectors.toList());

        final List<SomaticPeak> somaticPeaks = Lists.newArrayList();
        double initialWeight = positiveWeight(weightedVariants);

        double maxWeightPeak = 0;
        double maxWeightPeakWeight = 0;

        for(int i = 0; i < MAX_ITERATIONS; i++)
        {
            // Calculate peak
            double peak = mHistogram.peakPloidy(mPeakWidth, weightedVariants);

            if(peak <= 0)
                break;

            double offset = offset(peak);
            final WeightedPloidyHistogram peakHistogramFactory = new WeightedPloidyHistogram(mMaxVaf, mBinWidth, offset);
            final List<WeightedPloidy> peakPloidies = peakPloidies(peak, weightedVariants);
            double peakAverageWeight = averageWeight(peakPloidies);
            double[] peakHistogram = modelPeakHistogram(peak, peakPloidies);

            // sum and subtract modelled weight
            double[] currentHistogram = peakHistogramFactory.histogram(weightedVariants);
            double peakTotalWeight = Arrays.stream(peakHistogram).filter(x -> x > 0).sum();

            if(peakTotalWeight > maxWeightPeakWeight)
            {
                maxWeightPeakWeight = peakTotalWeight;
                maxWeightPeak = peak;
            }

            for(final ModifiableWeightedPloidy variant : weightedVariants)
            {
                int bucket = peakHistogramFactory.bucket(variant.ploidy());
                double currentWeight = variant.weight();
                double bucketWeight = currentHistogram[bucket];
                double peakWeight = peakHistogram[bucket];

                if(!Doubles.isZero(bucketWeight))
                {
                    double assignedWeight = abs(peakWeight / bucketWeight);
                    double newWeight = Doubles.isZero(bucketWeight) ? 0 : currentWeight - assignedWeight;
                    variant.setWeight(newWeight);
                }
            }

            if(peakTotalWeight < mMinPeakWeight || peakTotalWeight / initialWeight < mMinPeakPerc)
                break;

            somaticPeaks.add(ImmutableSomaticPeak.builder()
                    .alleleFrequency(peak)
                    .count((int)round(peakTotalWeight))
                    .build());

            // Decide if we should do another round
            double remainingWeight = positiveWeight(weightedVariants);
            double unexplainedWeight = remainingWeight / initialWeight;

            PPL_LOGGER.info(String.format("somatic peak(%.3f) weight(%.3f avg=%.3f) remaining(%.3f pct=%.3f)",
                    peak, peakTotalWeight, peakAverageWeight, remainingWeight, unexplainedWeight));

            if(Doubles.lessThan(unexplainedWeight, MAX_UNEXPLAINED_WEIGHT_PERCENT))
                break;
        }

        if(somaticPeaks.isEmpty())
        {
            PPL_LOGGER.info(String.format("using max somatic peak(%.3f weight=%.3f) not meeting criteria",
                    maxWeightPeak, maxWeightPeakWeight));
            return maxWeightPeak;
        }

        // return the highest VAF peak
        double maxVafPeak = somaticPeaks.stream()
                .filter(x -> x.alleleFrequency() <= 0.5)
                .mapToDouble(x -> x.alleleFrequency())
                .max().orElse(0);

        return maxVafPeak;
    }

    private double positiveWeight(@NotNull final List<? extends WeightedPloidy> weightedPloidies)
    {
        return weightedPloidies.stream().mapToDouble(x -> Math.max(0, x.weight())).sum();
    }

    double offset(double peak)
    {
        return peak - round(peak / mBinWidth) * mBinWidth;
    }

    @NotNull
    private List<WeightedPloidy> peakPloidies(double peak, @NotNull final List<? extends WeightedPloidy> allPloidies)
    {
        return allPloidies.stream()
                .filter(x -> Doubles.greaterThan(x.ploidy(), peak - mBinWidth / 2) && Doubles.lessThan(x.ploidy(),
                        peak + mBinWidth / 2))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    double[] modelPeakHistogram(double peak, @NotNull final List<WeightedPloidy> peakPloidies)
    {
        double offset = offset(peak);

        int maxBucket = bucket(mMaxVaf);
        double[] result = new double[maxBucket + 1];
        double[] weight = scalingFactor(peak, peakPloidies);

        int startBucket = bucket(peak - offset);

        // Forwards until unlikely...
        for(int i = startBucket; i <= maxBucket; i++)
        {
            double ploidy = i * mBinWidth + offset;
            double likelihood = likelihood(ploidy, weight, peakPloidies);
            result[i] = likelihood;
            if(Doubles.isZero(likelihood))
            {
                break;
            }
        }

        // Backwards until unlikely...
        for(int i = startBucket - 1; i >= 0; i--)
        {
            double ploidy = i * mBinWidth + offset;
            double likelihood = likelihood(ploidy, weight, peakPloidies);
            result[i] = likelihood;
            if(Doubles.isZero(likelihood))
            {
                break;
            }
        }

        return result;
    }

    private double likelihood(double ploidy, double[] scalingFactor, @NotNull List<WeightedPloidy> ploidies)
    {
        double result = 0;
        for(int i = 0; i < scalingFactor.length; i++)
        {
            result += scalingFactor[i] * ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    private double[] scalingFactor(double ploidy, @NotNull List<WeightedPloidy> ploidies)
    {
        double[] result = new double[ploidies.size()];
        for(int i = 0; i < ploidies.size(); i++)
        {
            result[i] = ploidies.get(i).weight() / ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    double ploidyLikelihood(double ploidy, @NotNull final WeightedPloidy weighted)
    {
        final String binomialKey = weighted.alleleReadCount() + ":" + weighted.totalReadCount();
        final BinomialDistribution binomialDistribution = mBinomialDistributionMap.computeIfAbsent(binomialKey,
                s -> new BinomialDistribution(weighted.totalReadCount(), weighted.alleleFrequency()));

        double lowerBoundAlleleReadCount = Math.max(0, ploidy - mBinWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int lowerBoundAlleleReadCountRounded = (int) round(lowerBoundAlleleReadCount);
        double lowerBoundAddition = lowerBoundAlleleReadCountRounded + 0.5 - lowerBoundAlleleReadCount;

        double upperBoundAlleleReadCount = Math.max(0, ploidy + mBinWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int upperBoundAlleleReadCountRounded = (int) round(upperBoundAlleleReadCount);
        double upperBoundSubtraction = upperBoundAlleleReadCountRounded + 0.5 - upperBoundAlleleReadCount;

        double rawResult =
                binomialDistribution.cumulativeProbability(upperBoundAlleleReadCountRounded) - binomialDistribution.cumulativeProbability(
                        lowerBoundAlleleReadCountRounded) + lowerBoundAddition * binomialDistribution.probability(
                        lowerBoundAlleleReadCountRounded) - upperBoundSubtraction * binomialDistribution.probability(
                        upperBoundAlleleReadCountRounded);

        return round(rawResult * 100) / 100d;
    }

    private int bucket(double ploidy)
    {
        return (int) round(ploidy / mBinWidth);
    }

    private static double averageWeight(@NotNull final List<WeightedPloidy> ploidies)
    {
        int count = ploidies.size();
        if(count == 0)
        {
            return 0;
        }

        return ploidies.stream().mapToDouble(WeightedPloidy::weight).sum() / count;
    }

    private static final double UPPER_PEAK_PROB = 0.95;

    public static double calcProbabilityUpperBound(int varCount, double peak)
    {
        final BinomialDistribution binomialDistribution = new BinomialDistribution(varCount, peak);

        // find the level above the peak with a probability of 95%
        int iterations = 0;

        int expCount = (int)round(peak * varCount);
        int lowerBound = expCount;
        int upperBound = (int)round(0.5 * varCount);
        int currentCount = (int)round((lowerBound + upperBound) * 0.5);

        while(iterations < 20)
        {
            double prob = binomialDistribution.cumulativeProbability(currentCount);

            if(abs(prob - UPPER_PEAK_PROB) < 0.005)
                break;

            if(prob < UPPER_PEAK_PROB)
            {
                // raise the count estimate
                if(currentCount >= upperBound - 1)
                    break;

                lowerBound = currentCount;
                currentCount = (int)round((currentCount + upperBound) * 0.5);
            }
            else
            {
                // lower the count
                if(currentCount <= lowerBound + 1)
                    break;

                upperBound = currentCount;
                currentCount = (int)round((currentCount + lowerBound) * 0.5);
            }

            ++iterations;
        }

        return currentCount / (double)varCount;
    }
}
