package com.hartwig.hmftools.common.variant.repeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RepeatContextFactoryTest
{
    @Test
    public void canDetermineRepeatContext()
    {
        String refGenome = "GATCGATCGATCGGAAAA";

        assertRepeatContext(3, "GATC", 5, refGenome);
        assertRepeatContext(3, "ATCG", 12, refGenome);
        assertRepeatContext(4, "A", 15, refGenome);
    }

    @Test
    public void testAdditionalBasesAtEnd()
    {
        String refGenome = "ATGTTTGTTTGTTTGAA";
        assertRepeatContext(1, 14, 3, "TGTT", 2, refGenome);
    }

    @Test
    public void testForwardsAndBackwardsCount()
    {
        assertRepeatContext(1, 6, 3, "CA", 5, "TCACACATTT");

        assertRepeatContext(0, 5, 3, "AC", 0, "ACACACTTTT");
        assertRepeatContext(0, 5, 3, "AC", 1, "ACACACTTTT");
        assertRepeatContext(0, 5, 3, "AC", 2, "ACACACTTTT");
        assertRepeatContext(0, 5, 3, "AC", 3, "ACACACTTTT");
        assertRepeatContext(0, 5, 3, "AC", 4, "ACACACTTTT");
        assertRepeatContext(0, 5, 3, "AC", 5, "ACACACTTTT");
        assertRepeatContext(6, 9, 4, "T", 6, "ACACACTTTT");
        assertRepeatContext(6, 9, 4, "T", 7, "ACACACTTTT");
        assertRepeatContext(6, 9, 4, "T", 8, "ACACACTTTT");
        assertRepeatContext(6, 9, 4, "T", 9, "ACACACTTTT");
    }

    @Test
    public void snvRepeat()
    {
        String refGenome = "GATCCCCCCT";
        assertNone(0, refGenome);
        assertNone(1, refGenome);
        assertNone(2, refGenome);
        assertRepeatContext(3, 8, 6, "C", 3, refGenome);
        assertRepeatContext(3, 8, 6, "C", 4, refGenome);
        assertRepeatContext(3, 8, 6, "C", 5, refGenome);
        assertRepeatContext(3, 8, 6, "C", 6, refGenome);
        assertRepeatContext(3, 8, 6, "C", 7, refGenome);
        assertRepeatContext(3, 8, 6, "C", 8, refGenome);

        refGenome = "GGGGGGATC";
        assertRepeatContext(0, 5, 6, "G", 0, refGenome);
        assertRepeatContext(0, 5, 6, "G", 1, refGenome);
        assertRepeatContext(0, 5, 6, "G", 2, refGenome);
        assertRepeatContext(0, 5, 6, "G", 3, refGenome);
        assertRepeatContext(0, 5, 6, "G", 4, refGenome);
        assertRepeatContext(0, 5, 6, "G", 5, refGenome);
        assertRepeatContext(0, 5, 6, "G", 5, refGenome);
        assertNone(6, refGenome);
    }

    @Test
    public void testMatchRepeatInBytes()
    {
        final String sequence = "AAGATCATC";
        assertFalse(RepeatContextFactory.match(3, 3, 0, sequence.getBytes()));
        assertFalse(RepeatContextFactory.match(3, 3, 1, sequence.getBytes()));
        assertFalse(RepeatContextFactory.match(3, 3, 2, sequence.getBytes()));
        assertTrue(RepeatContextFactory.match(3, 3, 3, sequence.getBytes()));
        assertFalse(RepeatContextFactory.match(3, 3, 4, sequence.getBytes()));
        assertFalse(RepeatContextFactory.match(3, 3, 5, sequence.getBytes()));
        assertTrue(RepeatContextFactory.match(3, 3, 6, sequence.getBytes()));
        assertFalse(RepeatContextFactory.match(3, 3, 7, sequence.getBytes()));
    }

    @Test
    public void testForwardRepeats()
    {
        final String sequence = "AAGATCATC";
        assertEquals(1, RepeatContextFactory.forwardRepeats(6, 3, sequence.getBytes()));
        assertEquals(2, RepeatContextFactory.forwardRepeats(3, 3, sequence.getBytes()));
    }

    @Test
    public void testBackwardRepeats()
    {
        final String sequence = "AAGATCATCATC";
        assertEquals(1, RepeatContextFactory.backwardRepeats(6, 3, sequence.getBytes()));
        assertEquals(2, RepeatContextFactory.backwardRepeats(9, 3, sequence.getBytes()));
    }

    private static void assertRepeatContext(int expectedCount, @NotNull String expectedBases, int index, @NotNull String sequence)
    {
        Optional<RepeatContext> optRepeatContextGATC = RepeatContextFactory.repeats(index, sequence);
        RepeatContext repeatContextGATC = optRepeatContextGATC.get();
        assertEquals(expectedCount, repeatContextGATC.count());
        assertEquals(expectedBases, repeatContextGATC.sequence());
    }

    private static void assertRepeatContext(int start, int end, int count, @NotNull String expectedBases, int index,
            @NotNull String sequence)
    {
        Optional<RepeatContext> optRepeatContextGATC = RepeatContextFactory.repeats(index, sequence);
        RepeatContext repeatContextGATC = optRepeatContextGATC.get();
        assertEquals(count, repeatContextGATC.count());
        assertEquals(start, repeatContextGATC.startIndex());
        assertEquals(end, repeatContextGATC.endIndex());
        assertEquals(expectedBases, repeatContextGATC.sequence());
    }

    private static void assertNone(int index, @NotNull String sequence)
    {
        Optional<RepeatContext> optRepeatContextGATC = RepeatContextFactory.repeats(index, sequence);
        assertFalse(optRepeatContextGATC.isPresent());
    }
}
