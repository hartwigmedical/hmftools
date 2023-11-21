package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.common.TestUtils.createSageConfig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collection;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class BaseQualityRecalibrationTest
{
    @Test
    public void testBaseQualityCounts()
    {
        BqrRegionReader bqrCounter = new BqrRegionReader(
                createSageConfig(), null, null, new BaseQualityResults());

        bqrCounter.initialise(new ChrBaseRegion("1", 100, 300));

        int pos1 = 100;
        BqrKey key1 = createKey('A', 'G', 30, pos1);
        BqrKey key2 = createKey('A', 'A', 20, pos1);
        BqrKey key3 = createKey('A', 'G', 15, pos1); // a repeated alt

        BaseQualityData bqData1 = bqrCounter.getOrCreateBaseQualData(pos1, key1.Ref, key1.TrinucleotideContext);

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
        BqrKey key4 = createKey('C', 'G', 25, pos2); // another repeated alt
        BaseQualityData bqData2 = bqrCounter.getOrCreateBaseQualData(pos2, key4.Ref, key4.TrinucleotideContext);

        for(int i = 0; i < 4; ++i)
        {
            bqData2.processReadBase(key4.Alt, key4.Quality);
        }

        int pos3 = 200;
        BqrKey key5 = createKey('A', 'G', 20, pos3); // an alt but not repeated
        BaseQualityData bqData3 = bqrCounter.getOrCreateBaseQualData(pos3, key5.Ref, key5.TrinucleotideContext);

        bqData3.processReadBase(key5.Alt, key5.Quality);

        bqrCounter.run();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(key1)).findFirst().orElse(null);
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

    private BqrKey createKey(char ref, char alt, int qual, int pos)
    {
        byte[] context = new byte[] { 65,  (byte)ref, 65};
        return new BqrKey((byte)ref, (byte)alt, context, (byte)qual);
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
