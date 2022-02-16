package com.hartwig.hmftools.sage.quality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;

import org.junit.Test;

public class BaseQualityRecalibrationTest
{
    @Test
    public void testBaseQualityCounts()
    {
        BaseQualityRegionCounter bqrCounter = new BaseQualityRegionCounter(
                new SageConfig(), null, null, new ChrBaseRegion("1", 100, 1000));

        Map<Integer, Map<BaseQualityKey,Integer>> bqMap = bqrCounter.getQualityMap();

        Map<BaseQualityKey,Integer> loc1 = Maps.newHashMap();
        int pos1 = 100;
        BaseQualityKey key1 = createKey('A', 'G', 30, pos1);
        BaseQualityKey key2 = createKey('A', 'A', 20, pos1);
        BaseQualityKey key3 = createKey('A', 'G', 15, pos1); // a repeated alt
        loc1.put(key1, 1);
        loc1.put(key2, 10);
        loc1.put(key3, 3);
        bqMap.put(pos1, loc1);

        // repeated alt at different locations
        int pos2 = 150;
        Map<BaseQualityKey,Integer> loc2 = Maps.newHashMap();
        BaseQualityKey key4 = createKey('C', 'G', 25, pos2); // another repeated alt
        loc2.put(key4, 4);
        bqMap.put(pos2, loc2);

        int pos3 = 200;
        Map<BaseQualityKey,Integer> loc3 = Maps.newHashMap();
        BaseQualityKey key5 = createKey('A', 'G', 20, pos3); // an alt but not repeated
        loc3.put(key5, 1);

        bqMap.put(pos3, loc3);

        bqrCounter.produceRegionCounts();

        Collection<QualityCounter> qualityCounts = bqrCounter.getQualityCounts();

        QualityCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(key1)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key2)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(10, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key3)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key4)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key5)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(1, qc.count());
    }

    private BaseQualityKey createKey(char ref, char alt, int qual, int pos)
    {
        byte[] context = new byte[] { 65,  (byte)alt, 65};
        return new BaseQualityKey((byte)ref, (byte)alt, context, (byte)qual);
    }

    @Test
    public void testBaseQualityAdjustment()
    {
        assertEquals(40, BaseQualityRecalibration.recalibratedQual(9999, 1), 0.1);
        assertEquals(30, BaseQualityRecalibration.recalibratedQual(9990, 10), 0.1);
        assertEquals(20, BaseQualityRecalibration.recalibratedQual(9900, 100), 0.1);
        assertEquals(10, BaseQualityRecalibration.recalibratedQual(9000, 1000), 0.1);
    }
}
