package com.hartwig.hmftools.sage.quality;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class QualityRecalibrationFactoryTest
{

    @Test
    public void testQual()
    {
        assertEquals(40, QualityRecalibrationFactory.recalibratedQual(9999, 1), 0.1);
        assertEquals(30, QualityRecalibrationFactory.recalibratedQual(9990, 10), 0.1);
        assertEquals(20, QualityRecalibrationFactory.recalibratedQual(9900, 100), 0.1);
        assertEquals(10, QualityRecalibrationFactory.recalibratedQual(9000, 1000), 0.1);
    }
}
