package com.hartwig.hmftools.sage.quality;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

public class QualityRecalibrationMap
{
    // recalibration results per sampleId
    private final Map<QualityRecalibrationKey,QualityRecalibrationRecord> mMap;

    public QualityRecalibrationMap(final List<QualityRecalibrationRecord> records)
    {
        mMap = records.stream().collect(Collectors.toMap(x -> x.Key, x -> x));
    }

    public double quality(byte ref, byte alt, byte[] trinucleotideContext, byte qual)
    {
        final QualityRecalibrationKey key = new QualityRecalibrationKey(ref, alt, trinucleotideContext, qual);

        QualityRecalibrationRecord record = mMap.get(key);
        return record != null ? record.RecalibratedQuality : 1.0;
    }
}
