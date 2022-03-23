package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class QualityRecalibrationMap
{
    // recalibration results per sampleId
    private final Map<BaseQualityKey,QualityRecalibrationRecord> mMap;

    public QualityRecalibrationMap(final List<QualityRecalibrationRecord> records)
    {
        mMap = Maps.newHashMap();

        for(QualityRecalibrationRecord record : records)
        {
            if(mMap.containsKey(record.Key))
            {
                QualityRecalibrationRecord existing = mMap.get(record.Key);
                SG_LOGGER.error("duplicate key({}) with existing key({}) count({})", record.Key, existing.Key, existing.Count);
            }
            else
            {
                mMap.put(record.Key, record);
            }
        }

        // mMap = records.stream().collect(Collectors.toMap(x -> x.Key, x -> x));
    }

    public double quality(byte ref, byte alt, byte[] trinucleotideContext, byte qual)
    {
        final BaseQualityKey key = new BaseQualityKey(ref, alt, trinucleotideContext, qual);

        QualityRecalibrationRecord record = mMap.get(key);
        return record != null ? record.RecalibratedQuality : qual;
    }
}
