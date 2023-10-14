package com.hartwig.hmftools.common.utils;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;

import org.junit.Test;

public class IntegersTest
{
    @Test
    public void testMedianNoEntries()
    {
        assertEquals(0, Integers.medianPositiveValue(Lists.newArrayList()));
    }

    @Test
    public void testMedianOneEntry()
    {
        assertEquals(5, Integers.medianPositiveValue(Lists.newArrayList(5)));
    }

    @Test
    public void testMedianTwoEntry()
    {
        assertEquals(6, Integers.medianPositiveValue(Lists.newArrayList(7, 5)));
    }

}
