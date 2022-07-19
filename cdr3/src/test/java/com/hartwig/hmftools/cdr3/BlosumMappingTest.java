package com.hartwig.hmftools.cdr3;

import static junit.framework.TestCase.assertEquals;

import org.junit.Test;

public class BlosumMappingTest
{
    @Test
    public void testBlosumMapping()
    {
        BlosumMapping mapping = new BlosumMapping();
        assertEquals(-3, mapping.map('W', 'A'));
        assertEquals(-2, mapping.map('R', 'Y'));
    }
}