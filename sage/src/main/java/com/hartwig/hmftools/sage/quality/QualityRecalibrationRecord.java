package com.hartwig.hmftools.sage.quality;

public class QualityRecalibrationRecord
{
    public final QualityRecalibrationKey Key;

    public final int Count;
    public final double RecalibratedQuality;

    public QualityRecalibrationRecord(final QualityRecalibrationKey key, final int count, final double recalibratedQual)
    {
        Key = key;
        Count = count;
        RecalibratedQuality = recalibratedQual;
    }
}
