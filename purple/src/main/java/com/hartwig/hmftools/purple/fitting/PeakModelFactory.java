package com.hartwig.hmftools.purple.fitting;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PeakModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(PeakModelFactory.class);
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    private static final int PEAK_BIN_COUNT = 10;
    private static final double PEAK_BIN_WIDTH = 0.01;

    private static final double MIN_AVERAGE_WEIGHT = 0.4;
    private static final double CLONAL_PLOIDY = 0.85;
    private static final int MAX_ITERATIONS = 10;

    private static final double MAX_UNEXPLAINED_WEIGHT_PERCENT = 0.01;

    private final double maxPloidy;
    private final double modelWidth;
    private final WeightedPloidyHistogram preciseHistogramFactory;
    private final Map<String, BinomialDistribution> binomialDistributionMap;

    public PeakModelFactory(final double maxPloidy, final double modelWidth) {
        this.modelWidth = modelWidth;
        this.maxPloidy = maxPloidy;
        this.preciseHistogramFactory = new WeightedPloidyHistogram(maxPloidy, PEAK_BIN_WIDTH);
        this.binomialDistributionMap = Maps.newHashMap();
    }

    @NotNull
    public List<PeakModel> model(@NotNull final List<ModifiableWeightedPloidy> weightedPloidies) {
        boolean hasValidSubclonalPeaks = false;
        final WeightedPloidyHistogram residualHistogram = new WeightedPloidyHistogram(0.85, modelWidth);
        double[] residualHistogramActual = residualHistogram.histogram(weightedPloidies);

        final List<ModifiablePeakModel> peakModel = Lists.newArrayList();
        double initialWeight = positiveWeight(weightedPloidies);

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            // Calculate peak
            double peak = preciseHistogramFactory.peakPloidy(PEAK_BIN_COUNT, weightedPloidies);
            double offset = offset(peak);
            final WeightedPloidyHistogram peakHistogramFactory = new WeightedPloidyHistogram(maxPloidy, modelWidth, offset);
            final List<WeightedPloidy> peakPloidies = peakPloidies(peak, weightedPloidies);
            double peakAverageWeight = averageWeight(peakPloidies);
            double[] peakHistogram = modelPeakHistogram(peak, peakPloidies);

            // Subtract modelled weight
            double[] currentHistogram = peakHistogramFactory.histogram(weightedPloidies);
            for (final ModifiableWeightedPloidy ploidy : weightedPloidies) {
                int bucket = peakHistogramFactory.bucket(ploidy.ploidy());
                double currentWeight = ploidy.weight();
                double bucketWeight = currentHistogram[bucket];
                double peakWeight = peakHistogram[bucket];
                double newWeight = Doubles.isZero(bucketWeight) ? 0 : currentWeight - Math.abs(peakWeight / bucketWeight);
                ploidy.setWeight(newWeight);
            }

            // Add results
            boolean isValidPeak = Doubles.greaterOrEqual(peakAverageWeight, MIN_AVERAGE_WEIGHT) && Doubles.greaterThan(peak, 0);
            boolean isSubclonal = Doubles.lessThan(peak, CLONAL_PLOIDY);
            hasValidSubclonalPeaks |= (isSubclonal && isValidPeak);
            for (int bucket = 0; bucket < peakHistogram.length; bucket++) {
                final ModifiablePeakModel model = ModifiablePeakModel.create()
                        .setBucket(bucket * modelWidth)
                        .setPeak(peak)
                        .setBucketWeight(peakHistogram[bucket])
                        .setPeakAvgWeight(peakAverageWeight)
                        .setIsSubclonal(isSubclonal)
                        .setIsValid(isValidPeak);
                peakModel.add(model);
            }

            // Decide if we should do another round
            double remainingWeight = positiveWeight(weightedPloidies);
            double unexplainedWeight = remainingWeight / initialWeight;

            LOGGER.debug("Peak: {}, Offset: {}, PeakAvgWeight: {}, Unexplained: {}",
                    new Object[] { FORMAT.format(peak), FORMAT.format(offset), FORMAT.format(peakAverageWeight),
                            FORMAT.format(unexplainedWeight) });

            if (Doubles.lessThan(unexplainedWeight, MAX_UNEXPLAINED_WEIGHT_PERCENT)) {
                break;
            }
        }

        // Scale results
        double totalModelWeight = peakModel.stream().filter(PeakModel::isValid).mapToDouble(PeakModel::bucketWeight).sum();
        double weightScalingFactor = initialWeight / totalModelWeight;
        LOGGER.debug("Weight scaling factor {}", String.format("%.4f", weightScalingFactor));

        final List<PeakModel> all = peakModel.stream().map(x -> x.setBucketWeight(x.bucketWeight() * 1)).collect(Collectors.toList());
        if (hasValidSubclonalPeaks) {
            return all;
        }

        // Find residual
        final List<PeakModel> validOnly = all.stream().filter(PeakModel::isValid).collect(Collectors.toList());
        final double[] residualHistogramModel = residualHistogram.modelHistogram(validOnly);
        all.addAll(residuals(residualHistogramActual, residualHistogramModel));

        return all;
    }

    @NotNull
    private List<PeakModel> residuals(double[] residualHistogramActual, double[] residualHistogramModel) {
        List<PeakModel> result = Lists.newArrayList();

        for (int i = 0; i < residualHistogramActual.length; i++) {
            double actualWeight = residualHistogramActual[i];
            double modelWeight = residualHistogramModel[i];

            final double residualPercent;
            if (Doubles.isZero(actualWeight)) {
                residualPercent = 1;
            } else {
                residualPercent = (actualWeight - modelWeight) / actualWeight;
            }

            if (Doubles.greaterThan(residualPercent, 0)) {
                result.add(ModifiablePeakModel.create()
                        .setBucket(i * modelWidth)
                        .setPeak(0)
                        .setBucketWeight(residualPercent)
                        .setPeakAvgWeight(1)
                        .setIsSubclonal(true)
                        .setIsValid(true));
            }
        }

        return result;
    }

    private double positiveWeight(@NotNull final List<? extends WeightedPloidy> weightedPloidies) {
        return weightedPloidies.stream().mapToDouble(x -> Math.max(0, x.weight())).sum();
    }

    double offset(double peak) {
        return peak - Math.round(peak / modelWidth) * modelWidth;
    }

    @NotNull
    private List<WeightedPloidy> peakPloidies(double peak, @NotNull final List<? extends WeightedPloidy> allPloidies) {
        return allPloidies.stream()
                .filter(x -> Doubles.greaterThan(x.ploidy(), peak - modelWidth / 2) && Doubles.lessThan(x.ploidy(), peak + modelWidth / 2))
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    double[] modelPeakHistogram(double peak, @NotNull final List<WeightedPloidy> peakPloidies) {
        double offset = offset(peak);

        int maxBucket = bucket(maxPloidy);
        double[] result = new double[maxBucket + 1];
        double[] weight = scalingFactor(peak, peakPloidies);

        int startBucket = bucket(peak - offset);

        // Forwards until unlikely...
        for (int i = startBucket; i <= maxBucket; i++) {
            double ploidy = i * modelWidth + offset;
            double likelihood = likelihood(ploidy, weight, peakPloidies);
            result[i] = likelihood;
            if (Doubles.isZero(likelihood)) {
                break;
            }
        }

        // Backwards until unlikely...
        for (int i = startBucket - 1; i >= 0; i--) {
            double ploidy = i * modelWidth + offset;
            double likelihood = likelihood(ploidy, weight, peakPloidies);
            result[i] = likelihood;
            if (Doubles.isZero(likelihood)) {
                break;
            }
        }

        return result;
    }

    private double likelihood(double ploidy, double[] scalingFactor, @NotNull List<WeightedPloidy> ploidies) {
        double result = 0;
        for (int i = 0; i < scalingFactor.length; i++) {
            result += scalingFactor[i] * ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    private double[] scalingFactor(double ploidy, @NotNull List<WeightedPloidy> ploidies) {
        double[] result = new double[ploidies.size()];
        for (int i = 0; i < ploidies.size(); i++) {
            result[i] = ploidies.get(i).weight() / ploidyLikelihood(ploidy, ploidies.get(i));
        }

        return result;
    }

    double ploidyLikelihood(double ploidy, @NotNull final WeightedPloidy weighted) {
        final String binomialKey = weighted.alleleReadCount() + ":" + weighted.totalReadCount();
        final BinomialDistribution binomialDistribution = binomialDistributionMap.computeIfAbsent(binomialKey,
                s -> new BinomialDistribution(weighted.totalReadCount(), weighted.alleleFrequency()));

        double lowerBoundAlleleReadCount = Math.max(0, ploidy - modelWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int lowerBoundAlleleReadCountRounded = (int) Math.round(lowerBoundAlleleReadCount);
        double lowerBoundAddition = lowerBoundAlleleReadCountRounded + 0.5 - lowerBoundAlleleReadCount;

        double upperBoundAlleleReadCount = Math.max(0, ploidy + modelWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int upperBoundAlleleReadCountRounded = (int) Math.round(upperBoundAlleleReadCount);
        double upperBoundSubtraction = upperBoundAlleleReadCountRounded + 0.5 - upperBoundAlleleReadCount;

        double rawResult =
                binomialDistribution.cumulativeProbability(upperBoundAlleleReadCountRounded) - binomialDistribution.cumulativeProbability(
                        lowerBoundAlleleReadCountRounded) + lowerBoundAddition * binomialDistribution.probability(
                        lowerBoundAlleleReadCountRounded) - upperBoundSubtraction * binomialDistribution.probability(
                        upperBoundAlleleReadCountRounded);

        return Math.round(rawResult * 100) / 100d;
    }

    private int bucket(double ploidy) {
        return (int) Math.round(ploidy / modelWidth);
    }

    private static double averageWeight(@NotNull final List<WeightedPloidy> ploidies) {
        int count = ploidies.size();
        if (count == 0) {
            return 0;
        }

        return ploidies.stream().mapToDouble(WeightedPloidy::weight).sum() / count;
    }
}
