package com.hartwig.hmftools.common.purple.repeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RepeatContextTest {

    @Test
    public void testSymmetric() {
        assertRepeats("GATACA", 6, "GATACAGATACAGATACAGATACAGATACAGATACA");
    }

    @Test
    public void testNonSymmetric() {
        final String sequence = "AAAAAAAAACA" + "BBBB";
        assertRepeats("A", 9, RepeatContextFactory.repeats(0, sequence));
        assertRepeats("A", 9, RepeatContextFactory.repeats(8, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory.repeats(9, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory.repeats(10, sequence));
        assertRepeats("B", 4, RepeatContextFactory.repeats(11, sequence));
    }

    private static void assertRepeats(final String expectedSequence, final int expectedCount, @NotNull final String victim) {
        for (int i = 0; i < victim.length(); i++) {
            assertRepeats(expectedSequence, expectedCount, RepeatContextFactory.repeats(i, victim));
        }
    }

    private static void assertRepeats(final String expectedSequence, final int expectedCount, @NotNull final Optional<RepeatContext> victim) {
        assertTrue(victim.isPresent());
        assertEquals(expectedSequence, victim.get().sequence());
        assertEquals(expectedCount, victim.get().count());
    }
}
