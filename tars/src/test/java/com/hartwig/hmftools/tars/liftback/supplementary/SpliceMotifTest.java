package com.hartwig.hmftools.tars.liftback.supplementary;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SpliceMotifTest
{
    @Test
    public void testClassifiesMotifTiers()
    {
        // canonical GT-AG and its reverse-complement CT-AC (strand unknown at scan time)
        assertEquals(Tier.CANONICAL, SpliceMotif.classify(bases("GT"), bases("AG")));
        assertEquals(Tier.CANONICAL, SpliceMotif.classify(bases("CT"), bases("AC")));

        // semi-canonical GC-AG / AT-AC and their reverse-complements CT-GC / GT-AT
        assertEquals(Tier.SEMI_CANONICAL, SpliceMotif.classify(bases("GC"), bases("AG")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceMotif.classify(bases("CT"), bases("GC")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceMotif.classify(bases("AT"), bases("AC")));
        assertEquals(Tier.SEMI_CANONICAL, SpliceMotif.classify(bases("GT"), bases("AT")));

        // lowercase normalises to the same tier
        assertEquals(Tier.CANONICAL, SpliceMotif.classify(bases("gt"), bases("ag")));
        assertEquals(Tier.CANONICAL, SpliceMotif.classify(bases("ct"), bases("ac")));

        // non-motif flanks, including a matching donor with a non-matching acceptor
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("AA"), bases("GG")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("NN"), bases("NN")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("GT"), bases("CC")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("CC"), bases("AG")));
    }

    @Test
    public void testRejectsNullOrWrongLengthInputs()
    {
        assertEquals(Tier.NONE, SpliceMotif.classify(null, bases("AG")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("GT"), null));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("G"), bases("AG")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("GT"), bases("A")));
        assertEquals(Tier.NONE, SpliceMotif.classify(bases("GTC"), bases("AG")));
    }
}
