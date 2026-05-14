package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

// covers LiftedAlignment.cigarHasRealNJunction — the anchor-aware N detector used by the
// discriminator. A 1-3 bp M anchor on either side of an N is the failure mode that previously
// caused BOTH_TX_JUNCTION_REF_MATCH to swap off a clean ref full-match.
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
        // the chr9 KLF4 case: 128M before N is fine, but 2M after N is not a real junction
        assertFalse(alignmentWithCigar("128M102N2M21S").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("147M3309N1M3S").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("130M172N3M18S").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testTinyLeadingAnchorRejected()
    {
        assertFalse(alignmentWithCigar("2M102N128M21S").cigarHasRealNJunction(MIN_ANCHOR));
        // soft-clip in front of the leading M doesn't rescue a too-short M anchor
        assertFalse(alignmentWithCigar("20S5M102N128M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testThresholdBoundary()
    {
        // exactly 8M on both sides — accepted
        assertTrue(alignmentWithCigar("8M100N8M").cigarHasRealNJunction(MIN_ANCHOR));
        // 7M on one side — rejected
        assertFalse(alignmentWithCigar("8M100N7M").cigarHasRealNJunction(MIN_ANCHOR));
        assertFalse(alignmentWithCigar("7M100N8M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testMultipleNJunctionsAllMustPass()
    {
        // every N must have valid anchors on both sides
        assertTrue(alignmentWithCigar("20M100N20M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
        // second junction has 5M middle anchor — rejected
        assertFalse(alignmentWithCigar("20M100N5M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
        // first junction has 5M middle anchor — rejected
        assertFalse(alignmentWithCigar("5M100N5M100N50M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testNonAdjacentMNotAnchor()
    {
        // 10M then 5I then N — the I breaks adjacency, so the N's leading anchor is 0
        assertFalse(alignmentWithCigar("10M5I100N20M").cigarHasRealNJunction(MIN_ANCHOR));
        // D after N before the M — trailing anchor is 0
        assertFalse(alignmentWithCigar("20M100N5D20M").cigarHasRealNJunction(MIN_ANCHOR));
    }

    @Test
    public void testRealJunctionAccepted()
    {
        // clean 50M + intron + 50M
        assertTrue(alignmentWithCigar("50M1000N50M").cigarHasRealNJunction(MIN_ANCHOR));
        // realistic 75bp+30bp split read
        assertTrue(alignmentWithCigar("75M5000N76M").cigarHasRealNJunction(MIN_ANCHOR));
    }
}
