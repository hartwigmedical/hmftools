package com.hartwig.hmftools.common.qual;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.BASE_QUAL_MINIMUM;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BaseQualAdjustmentTest
{
    @Test
    public void testBaseQualAdjustment()
    {
        double rawQual = 37;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        rawQual = 25;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        rawQual = 11;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        assertEquals(BASE_QUAL_MINIMUM, BaseQualAdjustment.adjustBaseQual(8.99));

    }
}
