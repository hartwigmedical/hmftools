package com.hartwig.hmftools.sage.quality;

public class QualityScores
{
    public final double CalcBaseQuality;
    public final double RecalibratedBaseQuality;
    public final int ModifiedMapQuality;
    public final double ModifiedBaseQuality;
    public final double ModifiedQuality;

    public QualityScores(
            double calcBaseQuality, double recalibratedBaseQuality, int modifiedMapQuality,
            double modifiedBaseQuality, double modifiedQuality)
    {
        CalcBaseQuality = calcBaseQuality;
        RecalibratedBaseQuality = recalibratedBaseQuality;
        ModifiedMapQuality = modifiedMapQuality;
        ModifiedBaseQuality = modifiedBaseQuality;
        ModifiedQuality = modifiedQuality;
    }
}
