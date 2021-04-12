package com.hartwig.hmftools.sigs;

import static com.hartwig.hmftools.common.sigs.DataUtils.round;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testRoundFunction()
    {
        double value = 1234.56789;
        assertEquals(1234.6, round(value, 1), 0.00001);
        assertEquals(1234.568, round(value, 3), 0.00001);
        assertEquals(1230, round(value, -1), 0.00001);
        assertEquals(1000, round(value, -3), 0.00001);
    }

}
