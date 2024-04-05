package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.qual.BqrRecord;

public class BqrRecordMap
{
    // recalibration results per sampleId
    private final Map<BqrKey, BqrRecord> mMap;

    public BqrRecordMap(final List<BqrRecord> records)
    {
        mMap = Maps.newHashMap();

        for(BqrRecord record : records)
        {
            if(mMap.containsKey(record.Key))
            {
                BqrRecord existing = mMap.get(record.Key);
                SG_LOGGER.error("duplicate key({}) with existing key({}) count({})", record.Key, existing.Key, existing.Count);
            }
            else
            {
                mMap.put(record.Key, record);
            }
        }
    }

    public double getQualityAdjustment(byte ref, byte alt, byte[] trinucleotideContext, byte qual, BqrReadType readType)
    {
        final BqrKey key = new BqrKey(ref, alt, trinucleotideContext, qual, readType);

        BqrRecord record = mMap.get(key);
        return record != null ? record.RecalibratedQuality : qual;
    }
}
