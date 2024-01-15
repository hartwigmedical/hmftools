package com.hartwig.hmftools.common.qual;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BaseQualAdjustmentTest
{
    @Test
    public void testBaseQualAdjustment()
    {
        double rawQual = 37;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        // assertEquals(37, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        // assertEquals(36, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99), 0.01);

        rawQual = 25;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        // assertEquals(25, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        // assertEquals(24, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99), 0.01);

        rawQual = 11;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        // assertEquals(11, BaseQualAdjustment.adjustBaseQual(rawQual - 0.099));
        // assertEquals(10, BaseQualAdjustment.adjustBaseQual(rawQual - 1));
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.99), 0.01);

        assertEquals(0, BaseQualAdjustment.adjustBaseQual(8.99));

    }
}
