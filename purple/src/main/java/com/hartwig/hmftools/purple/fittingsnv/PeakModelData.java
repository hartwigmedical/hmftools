package com.hartwig.hmftools.purple.fittingsnv;

public class PeakModelData
{
    public final double Peak;
    public final double PeakAvgWeight;
    public final double BucketWeight;
    public final double Bucket;

    public final boolean IsValid;
    public final boolean IsSubclonal;

    public PeakModelData(
            final double peak, final double peakAvgWeight, final double bucket, final double bucketWeight,
            final boolean isValid, final boolean isSubclonal)
    {
        Peak = peak;
        Bucket = bucket;
        BucketWeight = bucketWeight;
        PeakAvgWeight = peakAvgWeight;
        IsValid = isValid;
        IsSubclonal = isSubclonal;
    }
}
