package com.hartwig.hmftools.redux.splice.rescue;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;

import org.junit.Test;

public class SpliceMotifTest
{
    private static byte[] b(final String s)
    {
        return s.getBytes(StandardCharsets.US_ASCII);
    }

    @Test
    public void testCanonicalGTagForwardStrand()
    {
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(b("GT"), b("AG")));
    }

    @Test
    public void testCanonicalCTacReverseStrand()
    {
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(b("CT"), b("AC")));
    }

    @Test
    public void testSemiCanonicalGCagForwardStrand()
    {
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(b("GC"), b("AG")));
    }

    @Test
    public void testSemiCanonicalCTgcReverseStrand()
    {
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(b("CT"), b("GC")));
    }

    @Test
    public void testSemiCanonicalATacForwardStrand()
    {
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(b("AT"), b("AC")));
    }

    @Test
    public void testSemiCanonicalGTatReverseStrand()
    {
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(b("GT"), b("AT")));
    }

    @Test
    public void testCaseInsensitive()
    {
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(b("gt"), b("ag")));
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(b("ct"), b("ac")));
    }

    @Test
    public void testNoMotif()
    {
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("AA"), b("GG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("NN"), b("NN")));
        // half-canonical (only donor matches) doesn't count
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("GT"), b("CC")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("CC"), b("AG")));
    }

    @Test
    public void testNullOrWrongLengthInputs()
    {
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(null, b("AG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("GT"), null));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("G"), b("AG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("GT"), b("A")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(b("GTC"), b("AG")));
    }
}
