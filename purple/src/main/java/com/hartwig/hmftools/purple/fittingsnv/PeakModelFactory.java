package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.MAX_UNEXPLAINED_WEIGHT_PERCENT;
import static com.hartwig.hmftools.purple.PurpleConstants.PEAK_BIN_CLONAL_PLOIDY;
import static com.hartwig.hmftools.purple.PurpleConstants.PEAK_BIN_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.PEAK_BIN_MIN_AVERAGE_WEIGHT;
import static com.hartwig.hmftools.purple.PurpleConstants.PEAK_BIN_WIDTH;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class PeakModelFactory
{
    private static final int MAX_ITERATIONS = 10;
    private static final double MAX_HISTOGRAM_PLOIDY = 0.85;

    private final double mMaxPloidy;
    private final double mModelWidth;
    private final WeightedPloidyHistogram mPreciseHistogramFactory;
    private final Map<String, BinomialDistribution> mBinomialDistributionMap;

    public PeakModelFactory(final double maxPloidy, final double modelWidth)
    {
        mModelWidth = modelWidth;
        mMaxPloidy = maxPloidy;
        mPreciseHistogramFactory = new WeightedPloidyHistogram(maxPloidy, PEAK_BIN_WIDTH);
        mBinomialDistributionMap = Maps.newHashMap();
    }

    public List<PeakModelData> model(final List<WeightedPloidy> weightedPloidies)
    {
        boolean hasValidSubclonalPeaks = false;
        final WeightedPloidyHistogram residualHistogram = new WeightedPloidyHistogram(MAX_HISTOGRAM_PLOIDY, mModelWidth);
        double[] residualHistogramActual = residualHistogram.histogram(weightedPloidies);

        final List<PeakModelData> peakModel = Lists.newArrayList();
        double initialWeight = positiveWeight(weightedPloidies);

        for(int i = 0; i < MAX_ITERATIONS; i++)
        {
            // Calculate peak
            double peak = mPreciseHistogramFactory.peakPloidy(PEAK_BIN_COUNT, weightedPloidies);
            double offset = offset(peak);
            final WeightedPloidyHistogram peakHistogramFactory = new WeightedPloidyHistogram(mMaxPloidy, mModelWidth, offset);
            final List<WeightedPloidy> peakPloidies = peakPloidies(peak, weightedPloidies);
            double peakAverageWeight = averageWeight(peakPloidies);
            double[] peakHistogram = modelPeakHistogram(peak, peakPloidies);

            // Subtract modelled weight
            double[] currentHistogram = peakHistogramFactory.histogram(weightedPloidies);

            for(final WeightedPloidy ploidy : weightedPloidies)
            {
                int bucket = peakHistogramFactory.bucket(ploidy.ploidy());
                double currentWeight = ploidy.weight();
                double bucketWeight = currentHistogram[bucket];
                double peakWeight = peakHistogram[bucket];
                double newWeight = Doubles.isZero(bucketWeight) ? 0 : currentWeight - Math.abs(peakWeight / bucketWeight);
                ploidy.Weight = newWeight;
            }

            // Add results
            boolean isValidPeak = Doubles.greaterOrEqual(peakAverageWeight, PEAK_BIN_MIN_AVERAGE_WEIGHT) && Doubles.greaterThan(peak, 0);
            boolean isSubclonal = Doubles.lessThan(peak, PEAK_BIN_CLONAL_PLOIDY);
            hasValidSubclonalPeaks |= (isSubclonal && isValidPeak);
            for(int bucket = 0; bucket < peakHistogram.length; bucket++)
            {
                peakModel.add(new PeakModelData(
                        peak, peakAverageWeight, bucket * mModelWidth, peakHistogram[bucket], isValidPeak, isSubclonal));
            }

            // Decide if we should do another round
            double remainingWeight = positiveWeight(weightedPloidies);
            double unexplainedWeight = remainingWeight / initialWeight;

            PPL_LOGGER.trace(format("peak(%.2f) offset(%.2f) peakAvgWeight(%.2f) unexplained(%.2f)",
                    peak, offset, peakAverageWeight, unexplainedWeight));

            if(Doubles.lessThan(unexplainedWeight, MAX_UNEXPLAINED_WEIGHT_PERCENT))
            {
                break;
            }
        }

        // Scale results
        double totalModelWeight = peakModel.stream().filter(x -> x.IsValid).mapToDouble(x -> x.BucketWeight).sum();
        double weightScalingFactor = initialWeight / totalModelWeight;

        PPL_LOGGER.trace(format("weight scaling factor %.4f", weightScalingFactor));

        if(hasValidSubclonalPeaks)
            return peakModel;

        // find residuals
        List<PeakModelData> all = peakModel.stream().collect(Collectors.toList());
        List<PeakModelData> validOnly = all.stream().filter(x -> x.IsValid).collect(Collectors.toList());
        final double[] residualHistogramModel = residualHistogram.modelHistogram(validOnly);
        all.addAll(residuals(residualHistogramActual, residualHistogramModel));
        return all;
    }

    private List<PeakModelData> residuals(double[] residualHistogramActual, double[] residualHistogramModel)
    {
        List<PeakModelData> result = Lists.newArrayList();

        for(int i = 0; i < residualHistogramActual.length; i++)
        {
            double actualWeight = residualHistogramActual[i];
            double modelWeight = residualHistogramModel[i];

            final double residualPercent;
            if(Doubles.isZero(actualWeight))
            {
                residualPercent = 1;
            }
            else
            {
                residualPercent = (actualWeight - modelWeight) / actualWeight;
            }

            if(Doubles.greaterThan(residualPercent, 0))
            {
                result.add(new PeakModelData(
                        0, 1, i * mModelWidth, residualPercent, true, true));
            }
        }

        return result;
    }

    private double positiveWeight(final List<? extends WeightedPloidy> weightedPloidies)
    {
        return weightedPloidies.stream().mapToDouble(x -> Math.max(0, x.Weight)).sum();
    }

    public double offset(double peak)
    {
        return peak - Math.round(peak / mModelWidth) * mModelWidth;
    }

    private List<WeightedPloidy> peakPloidies(double peak, final List<? extends WeightedPloidy> allPloidies)
    {
        return allPloidies.stream()
                .filter(x -> Doubles.greaterThan(x.Ploidy, peak - mModelWidth / 2) && Doubles.lessThan(x.Ploidy,
                        peak + mModelWidth / 2))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    public double[] modelPeakHistogram(double peak, final List<WeightedPloidy> peakPloidies)
    {
        double offset = offset(peak);

        int maxBucket = bucket(mMaxPloidy);
        double[] result = new double[maxBucket + 1];
        double[] weight = scalingFactor(peak, peakPloidies);

        int startBucket = bucket(peak - offset);

        // Forwards until unlikely...
        for(int i = startBucket; i <= maxBucket; i++)
        {
            double ploidy = i * mModelWidth + offset;
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
            double ploidy = i * mModelWidth + offset;
            double likelihood = likelihood(ploidy, weight, peakPloidies);
            result[i] = likelihood;
            if(Doubles.isZero(likelihood))
            {
                break;
            }
        }

        return result;
    }

    private double likelihood(double ploidy, double[] scalingFactor, List<WeightedPloidy> ploidies)
    {
        double result = 0;
        for(int i = 0; i < scalingFactor.length; i++)
        {
            result += scalingFactor[i] * ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    private double[] scalingFactor(double ploidy, List<WeightedPloidy> ploidies)
    {
        double[] result = new double[ploidies.size()];
        for(int i = 0; i < ploidies.size(); i++)
        {
            result[i] = ploidies.get(i).Weight / ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    public double ploidyLikelihood(double ploidy, final WeightedPloidy weighted)
    {
        final String binomialKey = weighted.AlleleReadCount + ":" + weighted.TotalReadCount;
        final BinomialDistribution binomialDistribution = mBinomialDistributionMap.computeIfAbsent(binomialKey,
                s -> new BinomialDistribution(weighted.TotalReadCount, weighted.alleleFrequency()));

        double lowerBoundAlleleReadCount = Math.max(0, ploidy - mModelWidth / 2d) / weighted.Ploidy * weighted.AlleleReadCount;
        int lowerBoundAlleleReadCountRounded = (int) Math.round(lowerBoundAlleleReadCount);
        double lowerBoundAddition = lowerBoundAlleleReadCountRounded + 0.5 - lowerBoundAlleleReadCount;

        double upperBoundAlleleReadCount = Math.max(0, ploidy + mModelWidth / 2d) / weighted.Ploidy * weighted.AlleleReadCount;
        int upperBoundAlleleReadCountRounded = (int) Math.round(upperBoundAlleleReadCount);
        double upperBoundSubtraction = upperBoundAlleleReadCountRounded + 0.5 - upperBoundAlleleReadCount;

        double rawResult =
                binomialDistribution.cumulativeProbability(upperBoundAlleleReadCountRounded) - binomialDistribution.cumulativeProbability(
                        lowerBoundAlleleReadCountRounded) + lowerBoundAddition * binomialDistribution.probability(
                        lowerBoundAlleleReadCountRounded) - upperBoundSubtraction * binomialDistribution.probability(
                        upperBoundAlleleReadCountRounded);

        return Math.round(rawResult * 100) / 100d;
    }

    private int bucket(double ploidy)
    {
        return (int) Math.round(ploidy / mModelWidth);
    }

    private static double averageWeight(final List<WeightedPloidy> ploidies)
    {
        int count = ploidies.size();
        if(count == 0)
            return 0;

        return ploidies.stream().mapToDouble(x -> x.Weight).sum() / count;
    }
}
