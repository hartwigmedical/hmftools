package com.hartwig.hmftools.common.genome.chromosome;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ContigComparatorTest
{
    @Test
    public void testComparison()
    {
        assertEquals(0, ContigComparator.INSTANCE.compare("1", "1"));
        assertEquals(0, ContigComparator.INSTANCE.compare("chr1", "1"));
        assertEquals(0, ContigComparator.INSTANCE.compare("1", "chr1"));
        assertEquals(0, ContigComparator.INSTANCE.compare("chr1", "chr1"));

        assertDifference("1", "2");
        assertDifference("22", "X");
        assertDifference("X", "Y");
        assertDifference("Y", "M");
    }

    private static void assertDifference(final String first, final String second)
    {
        assertEquals(-1, ContigComparator.INSTANCE.compare(first, second) > 0 ? 1 : -1);
        assertEquals(1, ContigComparator.INSTANCE.compare(second, first) > 0 ? 1 : -1);
    }
}
