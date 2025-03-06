package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChangeContextTest extends ReversePaveTestBase
{
    PaddedExon ec = new PaddedExon(1, "", "", "TTTAAACCCGGG", 100, "CATG", "TACG");
    ChangeContext cc1 = new ChangeContext(ec,3, 6, true, 10);
    ChangeContext cc2 = new ChangeContext(ec,3, 6, false, 10);

    @Test
    public void insertionPointTest()
    {
        assertEquals(102, cc1.insertionPoint());
        assertEquals(108, cc2.insertionPoint());
    }

    @Test
    public void refBasesTest()
    {
        assertEquals("GTTT", cc(ec, 0, 2, true).refBases());
        assertEquals("TTTA", cc(ec, 1, 3, true).refBases());
        assertEquals("TAAA", cc(ec, 3, 5, true).refBases());
        assertEquals("CGGG", cc(ec, 9, 11, true).refBases());

        assertEquals("GTTT", cc(ec, 9, 11, false).refBases());
        assertEquals("TAAA", cc(ec, 6, 8, false).refBases());
        assertEquals("CCGG", cc(ec, 1, 3, false).refBases());
        assertEquals("CGGG", cc(ec, 0, 2, false).refBases());
    }

    private ChangeContext cc(PaddedExon pe, int start, int end, boolean isPositiveStrand)
    {
        return new ChangeContext(pe, start, end, isPositiveStrand, 10);
    }
}
