package com.hartwig.hmftools.common.variant.repeat;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import org.junit.Test;

public class RepeatContextFactoryTest {

    @Test
    public void canDetermineRepeatContext() {
        String refGenome = "GATCGATCGATCGAAAAA";

        assertRepeatContext(3, "GATC", 5, refGenome);
        assertRepeatContext(3, "ATCG", 12, refGenome);
        assertRepeatContext(5, "A", 15, refGenome);
    }

    @Test
    public void snvRepeat() {
        String refGenome = "GATCCCCCC";
        assertEquals("", RepeatContextFactory.repeats(2, refGenome).map(RepeatContext::sequence).orElse(""));
        assertEquals("C", RepeatContextFactory.repeats(3, refGenome).map(RepeatContext::sequence).orElse(""));

        refGenome = "GGGGGGATC";
        assertEquals("G", RepeatContextFactory.repeats(5, refGenome).map(RepeatContext::sequence).orElse(""));
        assertEquals("", RepeatContextFactory.repeats(6, refGenome).map(RepeatContext::sequence).orElse(""));
    }

    @Test
    public void canCountForwardRepeats() {
        assertEquals(1, RepeatContextFactory.forwardRepeats("A", "A"));
        assertEquals(2, RepeatContextFactory.forwardRepeats("A", "AA"));
        assertEquals(3, RepeatContextFactory.forwardRepeats("A", "AAA"));

        assertEquals(0, RepeatContextFactory.forwardRepeats("AC", "A"));
        assertEquals(1, RepeatContextFactory.forwardRepeats("AC", "AC"));
        assertEquals(1, RepeatContextFactory.forwardRepeats("AC", "ACA"));
        assertEquals(2, RepeatContextFactory.forwardRepeats("AC", "ACAC"));
        assertEquals(2, RepeatContextFactory.forwardRepeats("AC", "ACACA"));
    }

    @Test
    public void canCountBackwardsRepeats() {
        assertEquals(1, RepeatContextFactory.backwardRepeats("A", "A"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("A", "AA"));
        assertEquals(3, RepeatContextFactory.backwardRepeats("A", "AAA"));

        assertEquals(0, RepeatContextFactory.backwardRepeats("AC", "A"));
        assertEquals(1, RepeatContextFactory.backwardRepeats("AC", "AC"));
        assertEquals(1, RepeatContextFactory.backwardRepeats("AC", "AAC"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("AC", "ACAC"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("AC", "AACAC"));
    }

    private static void assertRepeatContext(int expectedCount, String expectedBases, int index, String sequence) {
        Optional<RepeatContext> optRepeatContextGATC = RepeatContextFactory.repeats(index, sequence);
        RepeatContext repeatContextGATC = optRepeatContextGATC.get();
        assertEquals(expectedCount, repeatContextGATC.count());
        assertEquals(expectedBases, repeatContextGATC.sequence());

        optRepeatContextGATC = RepeatContextFactory.repeats(index, sequence.getBytes());
        repeatContextGATC = optRepeatContextGATC.get();
        assertEquals(expectedCount, repeatContextGATC.count());
        assertEquals(expectedBases, repeatContextGATC.sequence());
    }

}
