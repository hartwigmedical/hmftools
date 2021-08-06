package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class NoiseModel
{
    private final double mRequiredProbability;
    private final double mWeight;
    private final AminoAcidFrequency mAminoAcidFrequency;

    private final Map<Integer, Map<Integer,Double>> mValueCache;

    private final static int MAX_ITERATIONS = 10;
    private final static double MAX_PROB_DIFF_PERCENT = 0.1; // ie 10% difference

    public NoiseModel(final AminoAcidFrequency aminoAcidFrequency, final double probability, final double weight)
    {
        mRequiredProbability = probability;
        mWeight = weight;
        mAminoAcidFrequency = aminoAcidFrequency;
        mValueCache = Maps.newHashMap();
    }

    public boolean enabled() { return mRequiredProbability > 0 && mWeight < 1; }

    public int cacheSize()
    {
        return mValueCache.values().stream().mapToInt(x -> x.size()).sum();
    }

    public double getExpected(final char aminoAcid, int totalBinds, int observedBinds)
    {
        return getExpected(BindConstants.aminoAcidIndex(aminoAcid), totalBinds, observedBinds);
    }

    public double getExpected(final int aaIndex, int totalBinds, int observedBinds)
    {
        // calculate an expected frequency for a given probability but only if the observed count is less than expected
        // by the amino acid's expected frequency in the proteome
        double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aaIndex);
        double expectedFreq = aaFrequency * totalBinds;

        if(observedBinds >= expectedFreq)
            return observedBinds;

        int rndTotalBinds = roundCount(totalBinds);
        int rndObservedBinds = roundCount(observedBinds);

        Map<Integer,Double> cache = mValueCache.get(rndTotalBinds);

        if(cache != null)
        {
            Double cacheExpected = cache.get(rndObservedBinds);

            if(cacheExpected != null)
                return cacheExpected;
        }
        else
        {
            cache = Maps.newHashMap();
            mValueCache.put(rndTotalBinds, cache);
        }

        double expected = calcExpected(rndTotalBinds, rndObservedBinds, expectedFreq, mRequiredProbability);

        double weightedExpected = observedBinds * mWeight + (1 - mWeight) *  expected;

        cache.put(rndObservedBinds, weightedExpected);
        return weightedExpected;
    }

    public static double calcExpected(int totalBinds, int observedBinds, double expectedFreq, double reqProbability)
    {
        if(observedBinds >= expectedFreq)
            return observedBinds;

        double lowerValue = observedBinds;
        double upperValue = expectedFreq;
        double currentValue = (lowerValue + upperValue) * 0.5;

        int iterations = 0;

        while(iterations < MAX_ITERATIONS)
        {
            double currentRate = currentValue / totalBinds;
            BinomialDistribution distribution = new BinomialDistribution(totalBinds, currentRate);
            double prob = distribution.cumulativeProbability(observedBinds);

            double probDiff = abs(prob - reqProbability) / reqProbability;
            if(probDiff < MAX_PROB_DIFF_PERCENT)
                return currentValue;

            if(prob < reqProbability)
            {
                // expected is too high vs observed so lower it
                upperValue = currentValue;
            }
            else
            {
                lowerValue = currentValue;
            }

            double newValue = (lowerValue + upperValue) * 0.5;

            if(abs(newValue - currentValue) < 0.25)
                return currentValue;

            currentValue = newValue;
            ++iterations;
        }

        return currentValue;
    }

    private int roundCount(int totalBinds)
    {
        if(totalBinds < 100)
            return totalBinds;

        double tickSize;

        if(totalBinds < 1000)
            tickSize = 10;
        else if(totalBinds < 10000)
            tickSize = 100;
        else if(totalBinds < 100000)
            tickSize = 1000;
        else if(totalBinds < 1000000)
            tickSize = 10000;
        else
            tickSize = 100000;

        return (int)(round(totalBinds/tickSize) * tickSize);
    }

}
