package com.hartwig.hmftools.cobalt.calculations;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Assert;
import org.junit.Test;

public class GCPailsListTest
{
    @Test
    public void pailTest()
    {
        GCPailsList pailsList = new GCPailsList();
        assertEquals(new GCPail(0), pailsList.getGCPail(0.0001));
        assertEquals(new GCPail(0), pailsList.getGCPail(0.0002));
        assertEquals(new GCPail(0), pailsList.getGCPail(0.003));

        assertEquals(new GCPail(20), pailsList.getGCPail(0.2001));
        assertEquals(new GCPail(21), pailsList.getGCPail(0.21));
        assertEquals(new GCPail(57), pailsList.getGCPail(0.56999));
        assertEquals(new GCPail(99), pailsList.getGCPail(0.99499));
        assertEquals(new GCPail(100), pailsList.getGCPail(0.999));
        assertEquals(new GCPail(100), pailsList.getGCPail(1.0));
    }

    @Test
    public void getBucketsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        assertEquals(101, pailsList.getBuckets().size());
        for (int i = 0; i < 101; i++)
        {
            assertEquals(new GCPail(i), pailsList.getBuckets().get(i));
        }
    }
}
