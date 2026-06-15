package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

// covers LiftedAlignment.cigarHasRealNJunction — rejects micro-anchors (1-3 bp M on either side of N)
// that previously caused BOTH_TX_JUNCTION_REF_MATCH to swap off a clean ref full-match.
public class LiftedAlignmentTest
{
    private static final int MIN_ANCHOR = 8;

    private static LiftedAlignment alignmentWithCigar(final String cigar)
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "ctg", 1, cigar,
                "chr1", 100, cigar,
                0, 0,
                null, null, null,
                false, true);
    }

    @Test
    public void testNoNReturnsFalse()
    {
        assertFalse(alignmentWithCigar("151M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("100M51S").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testTinyTrailingAnchorRejected()
    {
        assertFalse(alignmentWithCigar("128M102N2M21S").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("147M3309N1M3S").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("130M172N3M18S").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testTinyLeadingAnchorRejected()
    {
        assertFalse(alignmentWithCigar("2M102N128M21S").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("20S5M102N128M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testThresholdBoundary()
    {
        assertTrue(alignmentWithCigar("8M100N8M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("8M100N7M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("7M100N8M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testMultipleNJunctionsAllMustPass()
    {
        assertTrue(alignmentWithCigar("20M100N20M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("20M100N5M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("5M100N5M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testNonAdjacentMNotAnchor()
    {
        // I/D between M and N breaks adjacency — anchor is 0
        assertFalse(alignmentWithCigar("10M5I100N20M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("20M100N5D20M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testRealJunctionAccepted()
    {
        assertTrue(alignmentWithCigar("50M1000N50M").cigarHasRealNJunction(MIN_ANCHOR));
        assertTrue(alignmentWithCigar("75M5000N76M").cigarHasRealNJunction(MIN_ANCHOR));
    }
}
