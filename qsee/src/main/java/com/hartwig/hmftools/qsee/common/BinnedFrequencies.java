package com.hartwig.hmftools.qsee.common;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.metrics.ValueFrequency;

public class BinnedFrequencies
{
    private final double[] mBinStarts;
    private final double[] mBinEnds;
    private final double[] mFrequencies;

    public BinnedFrequencies(double[] binStarts, double[] frequencies)
    {
        mBinStarts = binStarts;
        mBinEnds = createBinEndsRightOpen();
        mFrequencies = frequencies;
    }

    private double[] createBinEndsRightOpen()
    {
        double[] binEnds = new double[mBinStarts.length];

        System.arraycopy(mBinStarts, 1, binEnds, 0, mBinStarts.length-1);
        binEnds[mBinStarts.length-1] = Double.POSITIVE_INFINITY;

        return binEnds;
    }

    public static BinnedFrequencies fromValueFrequencies(List<ValueFrequency> valueFrequencies)
    {
        double[] binStarts = new double[valueFrequencies.size()];
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
    public double[] binStarts() { return mBinStarts; }
    public double[] binEnds() { return mBinEnds; }

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
}
