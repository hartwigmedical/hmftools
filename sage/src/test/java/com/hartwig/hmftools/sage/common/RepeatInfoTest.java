package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.RepeatBoundaries.findRepeatBoundaries;
import static com.hartwig.hmftools.sage.common.RepeatInfo.addIfUnique;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMaxRepeat;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMultiBaseRepeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class RepeatInfoTest
{
    @Test
    public void testRepeats()
    {
        //              0123456789
        String bases = "AAACCTTTTT";

        // first check limits
        RepeatInfo repeatInfo = findMultiBaseRepeat(bases.getBytes(), 6, 1, 4);
        assertNotNull(repeatInfo);
        assertEquals("T", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 7, 1, 4);
        assertNull(repeatInfo);

        //       01234567890
        bases = "AACGGTACGGT";

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 1, 5, 2);
        assertNotNull(repeatInfo);
        assertEquals("ACGGT", repeatInfo.Bases);
        assertEquals(2, repeatInfo.Count);
        assertEquals(10, repeatInfo.totalLength());
        assertEquals(5, repeatInfo.repeatLength());

        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 2, 5, 2);
        assertNull(repeatInfo);
    }

    @Test
    public void testMaxRepeats()
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

    @Test
    public void testFindRepeatBoundaries()
    {
        //                  0123456789012345678901234567890123456789
        String readBases = "ACGTAGCTTGCGCGCGCCACACACTACGTAGCT";

        // test 1: a transitioning repeat where the lower repeat is the max

        int readCoreStart = 13;
        int readCoreEnd = 18;

        RepeatBoundaries repeatBoundaries = findRepeatBoundaries(
                readBases.getBytes(), readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

        assertNotNull(repeatBoundaries);
        assertEquals("GC", repeatBoundaries.MaxRepeat.Bases);
        assertEquals(4, repeatBoundaries.MaxRepeat.Count);
        assertEquals(8, repeatBoundaries.LowerIndex);
        assertEquals(24, repeatBoundaries.UpperIndex);


        // test 2: now where the upper repeat is the max but the second is still used to set the lower index boundary

        //           0123456789012345678901234567890123456789
        readBases = "ACGTAGCTTGCGCGCCACACACACTACGTAGCT";

        readCoreStart = 11;
        readCoreEnd = 16;

        repeatBoundaries = findRepeatBoundaries(
                readBases.getBytes(), readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

        assertNotNull(repeatBoundaries);
        assertEquals("CA", repeatBoundaries.MaxRepeat.Bases);
        assertEquals(4, repeatBoundaries.MaxRepeat.Count);
        assertEquals(8, repeatBoundaries.LowerIndex);
        assertEquals(24, repeatBoundaries.UpperIndex);

        // test 3: with partial repeats at the end which need to be included
        //           0123456789012345678901234567890123456789
        readBases = "AACTTTCTTTCTTTCTTAGGGGCGGGGCGGGGCGGGGTT";

        readCoreStart = 15;
        readCoreEnd = 19;

        repeatBoundaries = findRepeatBoundaries(
                readBases.getBytes(), readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

        assertNotNull(repeatBoundaries);
        assertEquals("G", repeatBoundaries.MaxRepeat.Bases);
        assertEquals(4, repeatBoundaries.MaxRepeat.Count);
        assertEquals(1, repeatBoundaries.LowerIndex);
        assertEquals(37, repeatBoundaries.UpperIndex);

        // test 4: 3 repeats, the inner one the longest but requiring use of the outer 2

        //                     10        20        30
        //           0123456789012345678901234567890123456789
        readBases = "ACTGAAAAAAAACAAAAAAAAAAACTATATATATATTACT";

        readCoreStart = 11;
        readCoreEnd = 25;

        repeatBoundaries = findRepeatBoundaries(
                readBases.getBytes(), readCoreStart, readCoreEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT);

        assertNotNull(repeatBoundaries);
        assertEquals("A", repeatBoundaries.MaxRepeat.Bases);
        assertEquals(11, repeatBoundaries.MaxRepeat.Count);
        assertEquals(3, repeatBoundaries.AllRepeats.size());
        assertEquals(3, repeatBoundaries.LowerIndex);
        assertEquals(36, repeatBoundaries.UpperIndex);

    }

    @Test
    public void testUniqueRepeats()
    {
        List<RepeatInfo> repeats = Lists.newArrayList();

        RepeatInfo repeat = new RepeatInfo(0, "A", 7);
        RepeatInfo repeat2 = new RepeatInfo(0, "AAA", 2);
        assertTrue(repeat.isMultipleOf(repeat2));

        repeat2 = new RepeatInfo(0, "AA", 3);
        assertTrue(repeat.isMultipleOf(repeat2));

        repeat = new RepeatInfo(0, "ACT", 4);
        repeat2 = new RepeatInfo(0, "ACTACT", 2);
        assertTrue(repeat.isMultipleOf(repeat2));

        repeat = new RepeatInfo(0, "ACT", 4);
        repeat2 = new RepeatInfo(0, "ACTAGT", 2);
        assertFalse(repeat.isMultipleOf(repeat2));

        // too long
        repeat = new RepeatInfo(0, "ACT", 4);
        repeat2 = new RepeatInfo(0, "ACTACT", 3);
        assertFalse(repeat.isMultipleOf(repeat2));

        // alternative starts for the same repeat
        repeat = new RepeatInfo(0, "ACT", 4);
        repeat2 = new RepeatInfo(1, "CTA", 3);
        RepeatInfo repeat3 = new RepeatInfo(2, "TAC", 3);
        assertTrue(repeat.isAlternativeStart(repeat2));
        assertTrue(repeat.isAlternativeStart(repeat3));

        // not if longer or starts outside first repeat's bounds
        repeat2 = new RepeatInfo(1, "CTA", 5);
        assertFalse(repeat.isAlternativeStart(repeat2));

        repeat2 = new RepeatInfo(3, "CTA", 3);
        assertFalse(repeat.isAlternativeStart(repeat2));

        // or different bases
        repeat2 = new RepeatInfo(1, "CTG", 3);
        assertFalse(repeat.isAlternativeStart(repeat2));

        assertTrue(addIfUnique(repeats, new RepeatInfo(0, "ACT", 3)));
        assertTrue(addIfUnique(repeats, new RepeatInfo(0, "A", 9)));
        assertFalse(addIfUnique(repeats, new RepeatInfo(0, "AA", 2)));
        assertFalse(addIfUnique(repeats, new RepeatInfo(1, "AA", 2)));
        assertFalse(addIfUnique(repeats, new RepeatInfo(1, "CTA", 3)));
    }
}
