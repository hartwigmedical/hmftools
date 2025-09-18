package com.hartwig.hmftools.cobalt.calculations;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Test;

public class GCPailTest
{
    @Test
    public void medianTest()
    {
        GCPail pail = new GCPail(50);
        assertEquals(-1.0, pail.median(), 0.0001);

        pail.addReading(34.5);
        assertEquals(34.5, pail.median(), 0.0001);

        pail.addReading(34.5);
        pail.addReading(38.5);
        assertEquals(34.5, pail.median(), 0.0001);

        pail.addReading(38.5);
        assertEquals((34.5 + 38.5)/2.0, pail.median(), 0.0001);
    }

    @Test
    public void equalsTest()
    {
        GCPail pail1 = new GCPail(50);
        GCPail pail2 = new GCPail(50);
        assertEquals(pail1, pail2);
        GCPail pail3 = new GCPail(60);
        assertNotEquals(pail1, pail3);
    }

    @Test
    public void hashCodeTest()
    {
        GCPail pail1 = new GCPail(50);
        GCPail pail2 = new GCPail(50);
        assertEquals(pail1.hashCode(), pail2.hashCode());
    }
}
