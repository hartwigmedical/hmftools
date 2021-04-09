package com.hartwig.hmftools.sigs.sim;

public class SimSigFactors {

    final public int SigId; // must match the col index in the signatures file

    final public String Name;

    final public double SampleProbability; // % of samples with the signature

    // a weighting on the random time component - 0 to 1, 1 meaning full effect, anything less lowering the impact
    final public double TimeFactor;

    final public int MedianCount; // median count from this signatures for all samples which have it

    final public double RateFactor; // a multiplier for expression of this sig

    public static final int SIG_FACTOR_CSV_ITEM_COUNT = 6;

    public SimSigFactors(final String[] items)
    {
        int item = 0;
        SigId = Integer.parseInt(items[item++]);
        Name = items[item++];
        SampleProbability = Double.parseDouble(items[item++]);
        TimeFactor = Double.parseDouble(items[item++]);
        MedianCount = Integer.parseInt(items[item++]);
        RateFactor = Double.parseDouble(items[item++]);
    }

    public boolean isValid()
    {
        return SampleProbability >= 0 && SampleProbability <= 1
            && TimeFactor >= 0 && TimeFactor <= 1
            && RateFactor >= 0;
    }

    public String toString()
    {
        return String.format("%d - %s: sampleProb(%.3f) time(%.2f) medianCount(%d) rate(%.3f)",
                SigId, Name, SampleProbability, TimeFactor, MedianCount, RateFactor);
    }
}
