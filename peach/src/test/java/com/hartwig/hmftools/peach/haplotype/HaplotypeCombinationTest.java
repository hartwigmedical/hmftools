package com.hartwig.hmftools.peach.haplotype;

import org.junit.Test;

import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

public class HaplotypeCombinationTest
{
    @Test
    public void testEquals()
    {
        HaplotypeCombination combination = new HaplotypeCombination(Map.of("*1", 2));
        assertEquals(new HaplotypeCombination(Map.of("*1", 2)), combination);
        assertNotEquals(new HaplotypeCombination(Map.of("*1", 3)), combination);
        assertNotEquals(new HaplotypeCombination(Map.of("*1", 1, "*7", 1)), combination);
    }
}
