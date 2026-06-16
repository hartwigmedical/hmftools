package com.hartwig.hmftools.tars.liftback.rescue;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SpliceMotifTest
{
    @Test
    public void testClassifiesMotifTiers()
    {
        // canonical GT-AG and its reverse-complement CT-AC (strand unknown at scan time)
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(bases("GT"), bases("AG")));
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(bases("CT"), bases("AC")));

        // semi-canonical GC-AG / AT-AC and their reverse-complements CT-GC / GT-AT
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(bases("GC"), bases("AG")));
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(bases("CT"), bases("GC")));
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(bases("AT"), bases("AC")));
        assertEquals(SpliceMotif.TIER_SEMI_CANONICAL, SpliceMotif.classify(bases("GT"), bases("AT")));

        // lowercase normalises to the same tier
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(bases("gt"), bases("ag")));
        assertEquals(SpliceMotif.TIER_CANONICAL, SpliceMotif.classify(bases("ct"), bases("ac")));

        // non-motif flanks, including a matching donor with a non-matching acceptor
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("AA"), bases("GG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("NN"), bases("NN")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("GT"), bases("CC")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("CC"), bases("AG")));
    }

    @Test
    public void testRejectsNullOrWrongLengthInputs()
    {
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(null, bases("AG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("GT"), null));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("G"), bases("AG")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("GT"), bases("A")));
        assertEquals(SpliceMotif.TIER_NONE, SpliceMotif.classify(bases("GTC"), bases("AG")));
    }
}
