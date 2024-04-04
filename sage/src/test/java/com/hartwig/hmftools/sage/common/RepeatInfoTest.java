package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.common.RepeatInfo.findMaxRepeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class RepeatInfoTest
{
    @Test
    public void testRepeats()
    {
        //              01234567890123456789
        String bases = "ACGTACGTTTTTTGTACGT";

        int maxLength = 5;
        int minCount = 3;

        int indexEnd = bases.length() - 1;

        RepeatInfo repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, false, -1);

        assertNotNull(repeat);
        assertEquals("T", repeat.Bases);
        assertEquals(7, repeat.Index);
        assertEquals(6, repeat.Count);

        repeat = findMaxRepeat(bases.getBytes(), 10, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("T", repeat.Bases);
        assertEquals(7, repeat.Index);
        assertEquals(6, repeat.Count);

        // longer repeats

        //       01234567890123456789
        bases = "ACGTACGAGAGAGGTACGT";

        repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GA", repeat.Bases);
        assertEquals(6, repeat.Index);
        assertEquals(3, repeat.Count);

        //       0123456789012345678901234567890123456789
        bases = "ACGTACACGGAAAGGAAAGGAAAGGAAAGGAAAGTACGT";

        repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GGAAA", repeat.Bases);
        assertEquals(8, repeat.Index);
        assertEquals(5, repeat.Count);

        // finding and extending the start
        repeat = findMaxRepeat(bases.getBytes(), 17, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GGAAA", repeat.Bases);
        assertEquals(8, repeat.Index);
        assertEquals(5, repeat.Count);

        // longest when there are lots to choose from

        //                 10        20        30        40        50
        //       012345678901234567890123456789012345678901234567890123456789
        bases = "AAAACGCGCGCGTTTTTGCTTTGCTTTGCTTTAACGTACGTACGTACGTCCDDCCCDDD";

        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 20);

        assertNotNull(repeat);
        assertEquals("TTTGC", repeat.Bases);
        assertEquals(14, repeat.Index);
        assertEquals(3, repeat.Count);

        // checks that the repeat crosses the required index
        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 40);

        assertNotNull(repeat);
        assertEquals("ACGT", repeat.Bases);
        assertEquals(33, repeat.Index);
        assertEquals(4, repeat.Count);

        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 8);

        assertNotNull(repeat);
        assertEquals("CG", repeat.Bases);
        assertEquals(4, repeat.Index);
        assertEquals(4, repeat.Count);
    }
}
