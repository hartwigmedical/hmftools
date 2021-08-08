package com.hartwig.hmftools.neo.bind;

import static org.apache.commons.math3.util.FastMath.log;

public class CalcConstants
{
    public final double PeptideLengthWeight;
    public final double AlleleWeight;
    public final double WeightExponent;
    public final double MaxAffinity;
    public final double BindingAffinityLow;
    public final double BindingAffinityHigh;
    public final double NoiseProbability;
    public final double NoiseWeight;
    public final double GlobalWeight;
    public final boolean ApplyScaledCount;

    public CalcConstants(
            double peptideLengthWeight, double alleleWeight, double weightExponent,
            double noiseProb, double noiseWeight, double globalWeight,
            double maxAffinity, double bindingAffinityLow, double bindingAffinityHigh, boolean applyScaledCount)
    {
        PeptideLengthWeight = peptideLengthWeight;
        AlleleWeight = alleleWeight;
        WeightExponent = weightExponent;
        NoiseProbability = noiseProb;
        NoiseWeight = noiseWeight;
        GlobalWeight = globalWeight;
        MaxAffinity = maxAffinity;
        BindingAffinityLow = bindingAffinityLow;
        BindingAffinityHigh = bindingAffinityHigh;
        ApplyScaledCount = applyScaledCount;
    }

    public double deriveLevelScore(final double affinity)
    {
        if(affinity >= MaxAffinity)
            return 0;

        if(affinity <= 0)
            return 1;

        return 1 - log(MaxAffinity, affinity);
    }

    public double deriveAffinityPercent(final double affinity)
    {
        if(affinity >= BindingAffinityHigh)
            return 0;

        if(affinity <= BindingAffinityLow)
            return 1;

        return (BindingAffinityHigh - affinity) / (BindingAffinityHigh - BindingAffinityLow);
    }
}
