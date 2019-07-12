package com.hartwig.hmftools.common.variant.clonality;

import java.util.Collection;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class WeightedPloidyHistogram {

    private final double maxPloidy;
    private final double binWidth;
    private final double offset;

    WeightedPloidyHistogram(final double maxPloidy, final double binWidth) {
        this(maxPloidy, binWidth, 0);
    }

    WeightedPloidyHistogram(final double maxPloidy, final double binWidth, final double offset) {
        this.maxPloidy = maxPloidy;
        this.binWidth = binWidth;
        this.offset = offset;
    }

    int bucket(double ploidy) {
        return (int) Math.round((ploidy - offset) / binWidth);
    }

    double ploidy(int bucket) {
        return bucket * binWidth + offset;
    }

    @NotNull
    public double[] histogram(@NotNull final Collection<? extends WeightedPloidy> ploidies) {
        int maxBucket = bucket(maxPloidy);
        double[] result = new double[maxBucket];

        for (final WeightedPloidy ploidy : ploidies) {

            int bucket = bucket(ploidy.ploidy());
            if (bucket <= maxBucket) {
                result[bucket] = result[bucket] + ploidy.weight();
            }
        }

        return result;
    }

    public double peakPloidy(int peakBinCount, @NotNull final Collection<? extends WeightedPloidy> ploidies) {
        return peakPloidy(peakBinCount, histogram(ploidies));
    }

    public double peakPloidy(int peakBinCount, @NotNull double[] histogram) {
        int peakBucket = peakBucket(peakBinCount, histogram);
        return ploidy(peakBucket);
    }

    public static int peakBucket(int bucketBuffer, @NotNull double[] histogram) {

        int maxBucket = -1;
        double maxBucketWeight = 0;

        for (int i = 0; i < histogram.length; i++) {
            if (Doubles.greaterThan(histogram[i], 0)) {

                double bucketWeight = 0;
                for (int j = Math.max(0, i - bucketBuffer); j <= Math.min(histogram.length - 1, i + bucketBuffer); j++) {
                    bucketWeight += histogram[j];
                }

                if (Doubles.greaterThan(bucketWeight, maxBucketWeight)) {
                    maxBucketWeight = bucketWeight;
                    maxBucket = i;
                }
            }
        }

        return maxBucket;
    }

}
