package com.hartwig.hmftools.neo.bind;

public class CalcConstants
{
    public final double PeptideLengthWeightFactor;
    public final double PeptideLengthWeightMax;

    public final double AlleleWeightFactor;
    public final double AlleleWeightMax;

    public final double NoiseProbability;
    public final double NoiseWeight;
    public final double GlobalWeight;


    public CalcConstants(
            double peptideLengthWeightFactor, double peptideLengthWeightMax, double alleleWeightFactor, double alleleWeightMax,
            double noiseProb, double noiseWeight, double globalWeight)
    {
        PeptideLengthWeightFactor = peptideLengthWeightFactor;
        PeptideLengthWeightMax = peptideLengthWeightMax;
        AlleleWeightFactor = alleleWeightFactor;
        AlleleWeightMax = alleleWeightMax;
        NoiseProbability = noiseProb;
        NoiseWeight = noiseWeight;
        GlobalWeight = globalWeight;
    }
}
