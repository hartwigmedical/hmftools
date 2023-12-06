package com.hartwig.hmftools.common.qual;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BaseQualAdjustmentTest
{
    @Test
    public void testBaseQualAdjustment()
    {
        double rawQual = 37;
        assertEquals(37, BaseQualAdjustment.adjustBaseQual(rawQual + 1));
        assertEquals(37, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        assertEquals(36, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(36, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99));

        rawQual = 25;
        assertEquals(25, BaseQualAdjustment.adjustBaseQual(rawQual + 1));
        assertEquals(25, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        assertEquals(24, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(24, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99));

        rawQual = 11;
        assertEquals(11, BaseQualAdjustment.adjustBaseQual(rawQual + 1));
        assertEquals(11, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        assertEquals(10, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(10, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99));

        assertEquals(0, BaseQualAdjustment.adjustBaseQual(8.99));

    }
}
