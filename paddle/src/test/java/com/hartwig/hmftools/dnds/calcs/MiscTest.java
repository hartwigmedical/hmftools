package com.hartwig.hmftools.dnds.calcs;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testExpectedDrivers()
    {
        DndsCv victim = new DndsCv(115, 3.4091126);
        assertEquals(72.78686, victim.expectedDrivers(103), 0.00001);
    }

}
