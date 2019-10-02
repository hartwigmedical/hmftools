package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReadContextTest {

    @Test
    public void testExactLengthOnBothSide() {
        final String expected = "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC";
        final ReadContext victim = new ReadContext(15, expected.getBytes(), 0, new byte[]{});
        assertTrue(victim.isComplete());
        assertEquals(expected, victim.toString());
    }

    @Test
    public void testAdditionalLengthOnBothSide() {
        final String expected = "ATG" + "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC" + "ATG";
        final ReadContext victim = new ReadContext(18, expected.getBytes(), 0, new byte[]{});
        assertTrue(victim.isComplete());
        assertEquals("AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC", victim.toString());
    }

    @Test
    public void testShortLeft() {
        final String bytes = "AATCCCCCCCCCCCCCCC";
        final ReadContext victim = new ReadContext(2, bytes.getBytes(), 0, new byte[]{});
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testShortRight() {
        final String bytes = "AAAAAAAAAAAAAAATCCC";
        final ReadContext victim = new ReadContext(15, bytes.getBytes(), 0, new byte[]{});
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testIncompleteStillMatches() {
        final ReadContext full = new ReadContext(15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes(), 0, new byte[]{});
        final ReadContext shortLeft = new ReadContext(2, "AATCCCCCCCCCCCCCCC".getBytes(), 0, new byte[]{});
        final ReadContext shortRight = new ReadContext(15, "AAAAAAAAAAAAAAATCCC".getBytes(), 0, new byte[]{});

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
        final ReadContext full = new ReadContext(15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes(), 0, new byte[]{});
        final ReadContext diffAtAlt = new ReadContext(15, "AAAAAAAAAAAAAAAGCCCCCCCCCCCCCCC".getBytes(), 0, new byte[]{});
        final ReadContext diffLeft = new ReadContext(15, "TAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes(), 0, new byte[]{});
        final ReadContext diffRight = new ReadContext(15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCG".getBytes(), 0, new byte[]{});
        final ReadContext diffBoth = new ReadContext(15, "TAAAAAAAAAAAAAATCCCCCCCCCCCCCCG".getBytes(), 0, new byte[]{});

        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffAtAlt));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffLeft));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffRight));
        assertEquals(ReadContext.ReadContextMatch.NONE,full.match(diffBoth));
    }


}
