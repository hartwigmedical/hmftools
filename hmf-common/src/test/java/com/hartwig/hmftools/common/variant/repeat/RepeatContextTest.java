package com.hartwig.hmftools.common.variant.repeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RepeatContextTest
{
    @Test
    public void testSymmetric()
    {
        assertRepeats("GATACA", 6, "GATACAGATACAGATACAGATACAGATACAGATACA");
    }

    @Test
    public void testNonSymmetric()
    {
        final String sequence = "AAAAAAAAACA" + "BBBB";
        assertRepeats("A", 9, RepeatContextFactory_.repeats(0, sequence));
        assertRepeats("A", 9, RepeatContextFactory_.repeats(8, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory_.repeats(9, sequence));
        assertEquals(Optional.empty(), RepeatContextFactory_.repeats(10, sequence));
        assertRepeats("B", 4, RepeatContextFactory_.repeats(11, sequence));
        assertRepeats("B", 4, RepeatContextFactory_.repeats(12, sequence));
        assertRepeats("B", 4, RepeatContextFactory_.repeats(13, sequence));
        assertRepeats("B", 4, RepeatContextFactory_.repeats(14, sequence));
    }

    private static void assertRepeats(final String expectedSequence, final int expectedCount, @NotNull final String victim)
    {
        for(int i = 0; i < victim.length(); i++)
        {
            assertRepeats(expectedSequence, expectedCount, RepeatContextFactory_.repeats(i, victim));
        }
    }

    private static void assertRepeats(final String expectedSequence, final int expectedCount,
            @NotNull final Optional<RepeatContext> victim)
    {
        assertTrue(victim.isPresent());
        assertEquals(expectedSequence, victim.get().sequence());
        assertEquals(expectedCount, victim.get().count());
    }
}
