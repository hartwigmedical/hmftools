package com.hartwig.hmftools.wisp.purity.variant;

public class PurityCalcData
{
    public double RawPurityEstimate;
    public double PurityEstimate;
    public double PurityRangeLow;
    public double PurityRangeHigh;
    public double LodPurityEstimate;
    public double DualProbability;

    public ClonalityData Clonality;

    /*
    public final double PurityProbability;
    */

    public PurityCalcData()
    {
        RawPurityEstimate = 0;
        PurityEstimate = 0;
        PurityRangeLow = 0;
        PurityRangeHigh = 0;
        LodPurityEstimate = 0;
        DualProbability = 0;

        Clonality = ClonalityData.NO_RESULT;

        /*
        PurityProbability = purityProbability;
        PurityRangeLow = purityRangeLow;
        PurityRangeHigh = purityRangeHigh;
        */
    }
}
