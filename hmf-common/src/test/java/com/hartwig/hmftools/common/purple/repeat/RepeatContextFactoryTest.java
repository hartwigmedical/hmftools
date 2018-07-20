package com.hartwig.hmftools.common.purple.repeat;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RepeatContextFactoryTest {

    @Test
    public void testForwardsRepeat() {
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
    public void testBackwardsRepeat() {
        assertEquals(1, RepeatContextFactory.backwardRepeats("A", "A"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("A", "AA"));
        assertEquals(3, RepeatContextFactory.backwardRepeats("A", "AAA"));

        assertEquals(0, RepeatContextFactory.backwardRepeats("AC", "A"));
        assertEquals(1, RepeatContextFactory.backwardRepeats("AC", "AC"));
        assertEquals(1, RepeatContextFactory.backwardRepeats("AC", "AAC"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("AC", "ACAC"));
        assertEquals(2, RepeatContextFactory.backwardRepeats("AC", "AACAC"));
    }
}
