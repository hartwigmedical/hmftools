package com.hartwig.hmftools.common.purple.repeat;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class RepeatContextTest {


    @Test
    public void testSymmetric() {
        assertRepeats("GATACA", 6, "GATACAGATACAGATACAGATACAGATACAGATACA");
    }

    @Test
    public void testNonSymmetric() {
        final String sequence = "AAAAAAAAABA" + "BBBB";
        assertRepeats("A", 9, RepeatContextFactory.repeats(0, sequence));
        assertRepeats("A", 9, RepeatContextFactory.repeats(8, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory.repeats(9, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory.repeats(10, sequence));
        assertRepeats("B", 4, RepeatContextFactory.repeats(11, sequence));

    }

    private void assertRepeats(final String expectedSequence, final int expectedCount, @NotNull final String victim) {
        for (int i = 0; i < victim.length(); i++) {
            assertRepeats(expectedSequence, expectedCount, RepeatContextFactory.repeats(i, victim));
        }
    }

    private void assertRepeats(final String expectedSequence, final int expectedCount, @NotNull final Optional<RepeatContext> victim) {
        assertEquals(expectedSequence, victim.get().sequence());
        assertEquals(expectedCount, victim.get().count());
    }

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
