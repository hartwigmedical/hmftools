package com.hartwig.hmftools.sage.count;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReadContextTest {

    @Test
    public void testExactLengthOnBothSide() {
        final String expected = "AAAAAAAAAATCCCCCCCCCC";
        final ReadContext victim = new ReadContext(10, expected.getBytes());
        assertTrue(victim.isComplete());
        assertEquals(expected, victim.toString());
    }

    @Test
    public void testAdditionalLengthOnBothSide() {
        final String expected = "ATG" + "AAAAAAAAAATCCCCCCCCCC" + "ATG";
        final ReadContext victim = new ReadContext(13, expected.getBytes());
        assertTrue(victim.isComplete());
        assertEquals("AAAAAAAAAATCCCCCCCCCC", victim.toString());
    }

    @Test
    public void testShortLeft() {
        final String bytes = "AATCCCCCCCCCC";
        final ReadContext victim = new ReadContext(2, bytes.getBytes());
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testShortRight() {
        final String bytes = "AAAAAAAAAATCCC";
        final ReadContext victim = new ReadContext(10, bytes.getBytes());
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testIncompleteStillMatches() {
        final ReadContext full = new ReadContext(10, "AAAAAAAAAATCCCCCCCCCC".getBytes());
        final ReadContext shortLeft = new ReadContext(2, "AATCCCCCCCCCC".getBytes());
        final ReadContext shortRight = new ReadContext(10, "AAAAAAAAAATCCC".getBytes());

        assertEquals(ReadContext.ReadContextMatch.FULL, full.match(full));

        assertEquals(ReadContext.ReadContextMatch.PARTIAL, full.match(shortLeft));
        assertEquals(ReadContext.ReadContextMatch.PARTIAL, full.match(shortRight));

        assertEquals(ReadContext.ReadContextMatch.PARTIAL,shortLeft.match(full));
        assertEquals(ReadContext.ReadContextMatch.PARTIAL,shortRight.match(full));

        assertEquals(ReadContext.ReadContextMatch.NONE,shortLeft.match(shortRight));
        assertEquals(ReadContext.ReadContextMatch.NONE,shortRight.match(shortLeft));
    }

    @Test
    public void testMatch() {
        final ReadContext full = new ReadContext(10, "AAAAAAAAAATCCCCCCCCCC".getBytes());
        final ReadContext diffAtAlt = new ReadContext(10, "AAAAAAAAAAGCCCCCCCCCC".getBytes());
        final ReadContext diffLeft = new ReadContext(10, "TAAAAAAAAATCCCCCCCCCC".getBytes());
        final ReadContext diffRight = new ReadContext(10, "AAAAAAAAAATCCCCCCCCCG".getBytes());
        final ReadContext diffBoth = new ReadContext(10, "AGAAAAAAAATCCCCCCCCGC".getBytes());

        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffAtAlt));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffLeft));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffRight));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffBoth));
    }


}
