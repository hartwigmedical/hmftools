package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.MAX_RECALIBRATED_BASE_QUAL;

import java.util.List;
import java.util.Map;
import static java.lang.Math.min;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.qual.BqrRecord;
import com.hartwig.hmftools.common.qual.BqrReadStrand;

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

    public double getQualityAdjustment(byte ref, byte alt, byte[] trinucleotideContext, byte qual, BqrReadType readType, BqrReadStrand readStrand)
    {
        final BqrKey key = new BqrKey(ref, alt, trinucleotideContext, qual, readType, readStrand);

        BqrRecord record = mMap.get(key);
        return record != null ? record.RecalibratedQuality : min(qual, MAX_RECALIBRATED_BASE_QUAL);
    }
}
