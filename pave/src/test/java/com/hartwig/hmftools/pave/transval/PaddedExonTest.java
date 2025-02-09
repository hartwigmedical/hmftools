package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PaddedExonTest extends TransvalTestBase
{
    PaddedExon ec = new PaddedExon("", "", "TTTAAACCCGGG", 100, "CATG");
    PaddedExon ec2 = new PaddedExon("A", "", "TTTAAACCCGG", 100, "CATG");
    PaddedExon ec3 = new PaddedExon("", "A", "TTAAACCCGGG", 100, "CATG");
    PaddedExon ec6 = new PaddedExon("TA", "AT", "TTAACCGG", 100, "CATG");

    @Test
    public void basesBetweenTest()
    {
        assertEquals("T", ec.basesBetween(0,0));
        assertEquals("TT", ec.basesBetween(0,1));
        assertEquals("TTTAAACCCGGG", ec.basesBetween(0,11));
    }

    @Test
    public void basesBefore()
    {
        assertEquals(0, ec.numberOfBasesInPreviousExon());
        assertEquals(1, ec2.numberOfBasesInPreviousExon());
        assertEquals(2, ec6.numberOfBasesInPreviousExon());
    }

    @Test
    public void codonLocationTest()
    {
        assertEquals(0, ec.codonLocationInExonBody(1));
        assertEquals(3, ec.codonLocationInExonBody(2));
        assertEquals(6, ec.codonLocationInExonBody(3));
        assertEquals(-1, ec2.codonLocationInExonBody(0));
        assertEquals(2, ec2.codonLocationInExonBody(1));
        assertEquals(5, ec2.codonLocationInExonBody(2));
        assertEquals(-2, ec6.codonLocationInExonBody(0));
        assertEquals(1, ec6.codonLocationInExonBody(1));
        assertEquals(4, ec6.codonLocationInExonBody(2));
    }

    @Test
    public void getSplitSequenceTest()
    {
        assertEquals(new SplitCodonSequence("GGG", null, 109), ec.getSplitSequenceForCodons(4,1));
        assertEquals(new SplitCodonSequence("CCC", null, 106), ec.getSplitSequenceForCodons(3,1));
        assertEquals(new SplitCodonSequence("CCCGGG", null, 106), ec.getSplitSequenceForCodons(3,2));
        assertEquals(new SplitCodonSequence("TTTAAA", null, 100), ec.getSplitSequenceForCodons(1,2));
        assertEquals(new SplitCodonSequence("A", "TTTAA", 100), ec2.getSplitSequenceForCodons(0,2));
        assertEquals(new SplitCodonSequence("TAAACC", null, 102), ec2.getSplitSequenceForCodons(1,2));
        assertEquals(new SplitCodonSequence("TA", "TTAACCG", 100), ec6.getSplitSequenceForCodons(0,3));
    }

    @Test
    public void splitAtEndTest()
    {
        assertEquals(new SplitCodonSequence("GG", "A", 109), ec3.getSplitSequenceForCodons(4,1));
        assertEquals(new SplitCodonSequence("CCGGG", "A", 106), ec3.getSplitSequenceForCodons(3,2));

    }

}
