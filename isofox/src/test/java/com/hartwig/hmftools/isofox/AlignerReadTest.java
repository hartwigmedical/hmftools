package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.FragmentAllocator.starFragmentCount;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class AlignerReadTest
{
    private static final double EPSILON = 1e-9;

    @Test
    public void testStarUniqueFullWeight()
    {
        assertEquals(1.0, starFragmentCount(255), EPSILON);
        assertEquals(1.0, starFragmentCount(4), EPSILON);
    }

    @Test
    public void testStarMultiMapTiers()
    {
        assertEquals(0.5, starFragmentCount(3), EPSILON);
        assertEquals(0.33, starFragmentCount(2), EPSILON);
        assertEquals(0.2, starFragmentCount(1), EPSILON);
        assertEquals(0.1, starFragmentCount(0), EPSILON);
    }
}
