package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChangeContextTest extends TransvalTestBase
{
    PaddedExon ec = new PaddedExon("", "", "TTTAAACCCGGG", 100, "CATG");
    PaddedExon ec2 = new PaddedExon("A", "", "TTTAAACCCGG", 100, "CATG");
    PaddedExon ec3 = new PaddedExon("", "A", "TTAAACCCGGG", 100, "CATG");
    PaddedExon ec6 = new PaddedExon("TA", "AT", "TTAACCGG", 100, "CATG");
    ChangeContext cc1 = new ChangeContext(ec,3, 6, true, 10);

    @Test
    public void applyDuplicationTest()
    {
        assertEquals(cr("FKKPG", "TTTAAAAAACCCGGG"), cc(ec, 3, 5, true).applyDuplication());
        assertEquals(cr("FKPAG", "TTTAAACCCGCCGGG"), cc(ec, 7, 9, true).applyDuplication());
    }

    private ChangeContext cc(PaddedExon pe, int start, int end, boolean isPositiveStrand)
    {
        return new ChangeContext(pe, start, end, isPositiveStrand, 10);
    }

    private ChangeResult cr(String aminoAcids, String bases)
    {
        return new ChangeResult(aaSeq(aminoAcids), bases);
    }
}
