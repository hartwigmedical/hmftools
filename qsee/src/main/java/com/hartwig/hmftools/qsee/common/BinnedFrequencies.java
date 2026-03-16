package com.hartwig.hmftools.qsee.common;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public class BinnedFrequencies
{
    private final long[] mBinStarts;
    private final long[] mBinEnds;
    private final double[] mFrequencies;

    public BinnedFrequencies(long[] binStarts, double[] frequencies)
    {
        mBinStarts = binStarts;
        mBinEnds = createBinEndsRightOpen();
        mFrequencies = frequencies;
    }

    private long[] createBinEndsRightOpen()
    {
        long[] binEnds = new long[mBinStarts.length];

        System.arraycopy(mBinStarts, 1, binEnds, 0, mBinStarts.length-1);
        binEnds[mBinStarts.length-1] = Long.MAX_VALUE;

        return binEnds;
    }

    public static BinnedFrequencies fromValueFrequencies(List<ValueFrequency> valueFrequencies)
    {
        long[] binStarts = new long[valueFrequencies.size()];
        double[] frequencies = new double[valueFrequencies.size()];

        for(int i = 0; i < valueFrequencies.size(); i++)
        {
            ValueFrequency valueFrequency = valueFrequencies.get(i);
            binStarts[i] = valueFrequency.Value;
            frequencies[i] = valueFrequency.Count;
        }

        return new BinnedFrequencies(binStarts, frequencies);
    }

    public double[] frequencies() { return mFrequencies; }
    public long[] binStarts() { return mBinStarts; }
    public long[] binEnds() { return mBinEnds; }

    public double[] calcFrequencyDensities()
    {
        double[] frequencyDensities = new double[mFrequencies.length];

        for(int i = 0; i < frequencyDensities.length; i++)
        {

            double frequency = mFrequencies[i];

            double frequencyDensity;

            boolean isLastBin = i == frequencyDensities.length-1;
            if(isLastBin)
            {
                frequencyDensities[i] = frequency;
                continue;
            }

            double binWidth = mBinEnds[i] - mBinStarts[i];
            frequencyDensity = frequency / binWidth;

            frequencyDensities[i] = frequencyDensity;
        }

        return frequencyDensities;
    }

    public double[] calcProportionalDensities()
    {
        double totalFrequency = Arrays.stream(mFrequencies).sum();

        double[] frequencyDensities = calcFrequencyDensities();
        double[] proportionalDensities = new double[mFrequencies.length];

        for(int i = 0; i < proportionalDensities.length; i++)
        {
            proportionalDensities[i] = frequencyDensities[i] / totalFrequency;
        }

        return proportionalDensities;
    }

    public List<Feature> formProportionalDensityFeatures(String binTitle, FeatureType featureType, SourceTool sourceTool)
    {
        double[] proportionalDensities = calcProportionalDensities();

        List<Feature> features = new ArrayList<>();

        for(int i = 0; i < proportionalDensities.length; i++)
        {
            String binString = String.valueOf(mBinStarts[i]);
            String featureName = MultiFieldStringBuilder.formSingleField(binTitle, binString);

            FeatureKey key = new FeatureKey(featureName, featureType, sourceTool);
            Feature feature = new Feature(key, proportionalDensities[i]);

            features.add(feature);
        }

        return features;
    }
}
