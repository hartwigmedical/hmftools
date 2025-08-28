package com.hartwig.hmftools.purple.fittingsnv;

import java.util.Collection;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class WeightedPloidyHistogram
{

    public final double MaxPloidy;
    public final double BinWidth;
    public final double Offset;

    public WeightedPloidyHistogram(final double maxPloidy, final double binWidth)
    {
        this(maxPloidy, binWidth, 0);
    }

    public WeightedPloidyHistogram(final double maxPloidy, final double binWidth, final double offset)
    {
        MaxPloidy = maxPloidy;
        BinWidth = binWidth;
        Offset = offset;
    }

    public int bucket(double ploidy)
    {
        return (int) Math.round((ploidy - Offset) / BinWidth);
    }

    public double ploidy(int bucket)
    {
        return bucket * BinWidth + Offset;
    }

    public double[] histogram(@NotNull final Collection<? extends WeightedPloidy> ploidies)
    {
        return histogram(ploidies, WeightedPloidy::ploidy, WeightedPloidy::weight);
    }

    public double[] modelHistogram(final List<PeakModelData> peakModels)
    {
        int maxBucket = bucket(MaxPloidy);
        double[] result = new double[maxBucket + 1];

        for(PeakModelData peakModelData : peakModels)
        {
            int bucket = bucket(peakModelData.Bucket);
            if(bucket <= maxBucket)
            {
                result[bucket] = result[bucket] + peakModelData.BucketWeight;
            }
        }

        return result;
    }

    private <T> double[] histogram(@NotNull final Collection<T> elements, Function<T, Double> ploidy, Function<T, Double> weight)
    {
        int maxBucket = bucket(MaxPloidy);
        double[] result = new double[maxBucket + 1];

        for(final T element : elements)
        {
            int bucket = bucket(ploidy.apply(element));
            if(bucket <= maxBucket)
            {
                result[bucket] = result[bucket] + weight.apply(element);
            }
        }

        return result;
    }

    public double peakPloidy(int peakBinCount, final Collection<? extends WeightedPloidy> ploidies)
    {
        return peakPloidy(peakBinCount, histogram(ploidies));
    }

    private double peakPloidy(int peakBinCount, final double[] histogram)
    {
        int peakBucket = peakBucket(peakBinCount, histogram);
        return ploidy(peakBucket);
    }

    public static int peakBucket(int bucketBuffer, final double[] histogram)
    {
        int maxBucket = -1;
        double maxBucketWeight = 0;

        for(int i = 0; i < histogram.length; i++)
        {
            if(Doubles.greaterThan(histogram[i], 0))
            {

                double bucketWeight = 0;

                for(int j = Math.max(0, i - bucketBuffer); j <= Math.min(histogram.length - 1, i + bucketBuffer); j++)
                {
                    bucketWeight += histogram[j];
                }

                if(Doubles.greaterThan(bucketWeight, maxBucketWeight))
                {
                    maxBucketWeight = bucketWeight;
                    maxBucket = i;
                }
            }
        }

        return maxBucket;
    }
}
