package com.hartwig.hmftools.neo.bind;

import static org.apache.commons.math3.util.FastMath.log;

public class CalcConstants
{
    public final double PeptideLengthWeight;
    public final double AlleleWeight;
    public final double WeightExponent;
    public final double NoiseProbability;
    public final double NoiseWeight;
    public final double GlobalWeight;

    public CalcConstants(
            double peptideLengthWeight, double alleleWeight, double weightExponent,
            double noiseProb, double noiseWeight, double globalWeight)
    {
        PeptideLengthWeight = peptideLengthWeight;
        AlleleWeight = alleleWeight;
        WeightExponent = weightExponent;
        NoiseProbability = noiseProb;
        NoiseWeight = noiseWeight;
        GlobalWeight = globalWeight;
    }
}
