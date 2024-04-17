package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.common.TestUtils.createSageConfig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.Collections;

import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;

import org.junit.Test;

public class BaseQualityRecalibrationTest
{
    private static final String SAMPLE_ID = "SAMPLE_ID";
    @Test
    public void testBaseQualityCounts()
    {
        SageConfig config = createSageConfig();

        BqrRegionReader bqrCounter = new BqrRegionReader(
                config, null, null, new BaseQualityResults(), new BqrRecordWriter(config, SAMPLE_ID));

        bqrCounter.initialise(new ChrBaseRegion("1", 100, 300), Collections.emptySet());

        int pos1 = 100;
        BqrKey key1 = createKey('A', 'G', 30, BqrReadType.NONE);
        BqrKey key2 = createKey('A', 'A', 20, BqrReadType.NONE);
        BqrKey key3 = createKey('A', 'G', 15, BqrReadType.NONE); // a repeated alt

        BaseQualityData bqData1 = bqrCounter.getOrCreateBaseQualData(pos1, key1.Ref, key1.TrinucleotideContext, BqrReadType.NONE);

        bqData1.processReadBase(key1.Alt, key1.Quality);

        for(int i = 0; i < 10; ++i)
        {
            bqData1.processReadBase(key2.Alt, key2.Quality);
        }

        for(int i = 0; i < 3; ++i)
        {
            bqData1.processReadBase(key3.Alt, key3.Quality);
        }

        // repeated alt at different locations
        int pos2 = 150;
        BqrKey key4 = createKey('C', 'G', 25, BqrReadType.NONE); // another repeated alt
        BaseQualityData bqData2 = bqrCounter.getOrCreateBaseQualData(pos2, key4.Ref, key4.TrinucleotideContext, BqrReadType.NONE);

        for(int i = 0; i < 4; ++i)
        {
            bqData2.processReadBase(key4.Alt, key4.Quality);
        }

        int pos3 = 200;
        BqrKey key5 = createKey('A', 'G', 20, BqrReadType.NONE); // an alt but not repeated
        BaseQualityData bqData3 = bqrCounter.getOrCreateBaseQualData(pos3, key5.Ref, key5.TrinucleotideContext, BqrReadType.NONE);

        bqData3.processReadBase(key5.Alt, key5.Quality);

        for(int i = 0; i < 9; ++i)
        {
            bqData3.processReadBase(key2.Ref, key2.Quality); // AF of 10% but count of 1 is permitted
        }

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(key1)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key2)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(19, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key3)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key4)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key5)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(1, qc.count());
    }

    @Test
    public void testBaseQualityCountsByReadType()
    {
        SageConfig config = createSageConfig();

        BqrRegionReader bqrCounter = new BqrRegionReader(
                config, null, null, new BaseQualityResults(), new BqrRecordWriter(config, SAMPLE_ID));

        bqrCounter.initialise(new ChrBaseRegion("1", 100, 300), Collections.emptySet());

        char refBase = 'G';
        char altBase = 'G';
        byte qual1 = 37;
        byte qual2 = 25;

        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.NONE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.NONE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.NONE, qual1);

        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.SINGLE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.SINGLE, qual1);

        addReadBaseQual(bqrCounter, 100, refBase, altBase, BqrReadType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 101, refBase, altBase, BqrReadType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 102, refBase, altBase, BqrReadType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 103, refBase, altBase, BqrReadType.DUAL, qual1);

        addReadBaseQual(bqrCounter, 101, refBase, altBase, BqrReadType.DUAL, qual2);
        addReadBaseQual(bqrCounter, 102, refBase, altBase, BqrReadType.DUAL, qual2);

        BqrKey keyNone = createKey(refBase, altBase, qual1, BqrReadType.NONE);
        BqrKey keySingle = createKey(refBase, altBase, qual1, BqrReadType.SINGLE);
        BqrKey keyDualQ1 = createKey(refBase, altBase, qual1, BqrReadType.DUAL);
        BqrKey keyDualQ2 = createKey(refBase, altBase, qual2, BqrReadType.DUAL);

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNone)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(3, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keySingle)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualQ1)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(4, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualQ2)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());
    }

    private static void addReadBaseQual(
            final BqrRegionReader bqrCounter, int position, char ref, char alt, final BqrReadType readType, byte quality)
    {
        byte[] context = new byte[] { 65,  (byte)ref, 65};
        bqrCounter.getOrCreateBaseQualData(position, (byte)ref, context, readType).processReadBase((byte)alt, quality);
    }

    private BqrKey createKey(char ref, char alt, int qual, final BqrReadType readType)
    {
        byte[] context = new byte[] { 65,  (byte)ref, 65};
        return new BqrKey((byte)ref, (byte)alt, context, (byte)qual, readType);
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
