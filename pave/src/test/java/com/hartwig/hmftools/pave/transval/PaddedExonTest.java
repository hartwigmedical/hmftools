package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PaddedExonTest extends TransvalTestBase
{
    PaddedExon ec = new PaddedExon(0,"", "", "TTTAAACCCGGG", 100, "CATG");
    PaddedExon ec2 = new PaddedExon(2, "A", "", "TTTAAACCCGG", 100, "CATG");
    PaddedExon ec3 = new PaddedExon(3, "", "A", "TTAAACCCGGG", 100, "CATG");
    PaddedExon ec6 = new PaddedExon(6, "TA", "AT", "TTAACCGG", 100, "CATG");

    @Test
    public void fromStrandCoordinatesTest()
    {
        assertEquals(0, ec.fromStrandCoordinates(100, true));
        assertEquals(10, ec.fromStrandCoordinates(110, true));
        assertEquals(5, ec2.fromStrandCoordinates(105, true));

        assertEquals(12, ec.fromStrandCoordinates(100, false));
        assertEquals(2, ec.fromStrandCoordinates(110, false));
        assertEquals(6, ec2.fromStrandCoordinates(105, false));
    }

    @Test
    public void toStrandCoordinatesTest()
    {
        assertEquals(100, ec.toStrandCoordinates(0, true));
        assertEquals(110, ec.toStrandCoordinates(10, true));
        assertEquals(105, ec2.toStrandCoordinates(5, true));

        assertEquals(112, ec.toStrandCoordinates(0, false));
        assertEquals(111, ec.toStrandCoordinates(1, false));
        assertEquals(102, ec.toStrandCoordinates(10, false));
        assertEquals(106, ec2.toStrandCoordinates(5, false));

        for(int i=0; i<12; i++)
        {
            assertEquals(i, ec.fromStrandCoordinates(ec.toStrandCoordinates(i, true),true));
            assertEquals(i, ec3.fromStrandCoordinates(ec3.toStrandCoordinates(i, true),true));
            assertEquals(i, ec.fromStrandCoordinates(ec.toStrandCoordinates(i, false),false));
            assertEquals(i, ec3.fromStrandCoordinates(ec3.toStrandCoordinates(i, false),false));
        }
    }

    @Test
    public void baseSequenceWithBasesReplacedTest()
    {
        assertEquals("GTGAAACCCGGG", ec.baseSequenceWithBasesReplaced(0, "GTG", true));
        assertEquals("TGTGAACCCGGG", ec.baseSequenceWithBasesReplaced(1, "GTG", true));
        assertEquals("TTTCCCCCCGGG", ec.baseSequenceWithBasesReplaced(3, "CCC", true));
        assertEquals("TTTAAACCCTAC", ec.baseSequenceWithBasesReplaced(9, "TAC", true));
        assertEquals("AGGGAAACCCGG", ec2.baseSequenceWithBasesReplaced(0, "GGG", true));
        assertEquals("TATTAACTTTAT", ec6.baseSequenceWithBasesReplaced(5, "TTT", true));

        assertEquals("GTAGGGTTTAAA", ec.baseSequenceWithBasesReplaced(0, "TAC", false));
        assertEquals("CCCGTATTTAAA", ec.baseSequenceWithBasesReplaced(3, "TAC", false));
        assertEquals("CCGAAATTAAAT", ec2.baseSequenceWithBasesReplaced(3, "TTT", false));
        assertEquals("ATCCGGTACGTA", ec6.baseSequenceWithBasesReplaced(5, "CGT", false));
    }

    @Test
    public void baseAtTest()
    {
        assertEquals("T", ec.baseAt(0, true));
        assertEquals("A", ec.baseAt(3, true));
        assertEquals("G", ec.baseAt(11, true));
        assertEquals("T", ec2.baseAt(0, true));
        assertEquals("T", ec2.baseAt(1, true));
        assertEquals("T", ec6.baseAt(0, true));
        assertEquals("A", ec6.baseAt(3, true));

        assertEquals("G", ec.baseAt(0, false));
        assertEquals("C", ec.baseAt(3, false));
        assertEquals("T", ec.baseAt(11, false));
        assertEquals("G", ec2.baseAt(0, false));
        assertEquals("G", ec2.baseAt(1, false));
        assertEquals("G", ec6.baseAt(0, false));
        assertEquals("C", ec6.baseAt(3, false));
    }

    @Test
    public void codonLocationTest()
    {
        assertEquals(0, ec.codonLocationInExonBody(1, true));
        assertEquals(3, ec.codonLocationInExonBody(2, true));
        assertEquals(6, ec.codonLocationInExonBody(3, true));
        assertEquals(-1, ec2.codonLocationInExonBody(0, true));
        assertEquals(2, ec2.codonLocationInExonBody(1, true));
        assertEquals(5, ec2.codonLocationInExonBody(2, true));
        assertEquals(-2, ec6.codonLocationInExonBody(0, true));
        assertEquals(1, ec6.codonLocationInExonBody(1, true));
        assertEquals(4, ec6.codonLocationInExonBody(2, true));

        assertEquals(9, ec.codonLocationInExonBody(1, false));
        assertEquals(6, ec.codonLocationInExonBody(2, false));
        assertEquals(3, ec.codonLocationInExonBody(3, false));
        assertEquals(11, ec2.codonLocationInExonBody(0, false));
        assertEquals(8, ec2.codonLocationInExonBody(1, false));
        assertEquals(5, ec2.codonLocationInExonBody(2, false));
        assertEquals(2, ec2.codonLocationInExonBody(3, false));
        assertEquals(9, ec3.codonLocationInExonBody(0, false));
        assertEquals(6, ec3.codonLocationInExonBody(1, false));
        assertEquals(3, ec3.codonLocationInExonBody(2, false));
        assertEquals(0, ec3.codonLocationInExonBody(3, false));
//        assertEquals(-2, ec6.codonLocationInExonBody(0, false));
//        assertEquals(1, ec6.codonLocationInExonBody(1, false));
//        assertEquals(4, ec6.codonLocationInExonBody(2, false));
    }

    @Test
    public void baseImmediatelyBeforeTest()
    {
        assertEquals("T", ec.baseImmediatelyBefore(1));
        assertEquals("T", ec.baseImmediatelyBefore(3));
        assertEquals("G", ec.baseImmediatelyBefore(11));

        assertEquals("G", ec2.baseImmediatelyBefore(0));
        assertEquals("T", ec2.baseImmediatelyBefore(1));

        assertEquals("G", ec6.baseImmediatelyBefore(0));
        assertEquals("T", ec6.baseImmediatelyBefore(1));
    }

    @Test
    public void getCodonTest()
    {
        CodonWithinExons cwe = ec.getCodon(2, true);
        assertEquals("AAA", cwe.codon());
        assertEquals( 103, cwe.strandLocationOfStartOfVariablePart());
        assertEquals("AAA", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix());

        cwe = ec2.getCodon(0, true);
        assertEquals("ATT", cwe.codon());
        assertEquals("TT", cwe.variablePart());
        assertEquals("A", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix());

        cwe = ec2.getCodon(1, true);
        assertEquals("TAA", cwe.codon());
        assertEquals("TAA", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix());

        cwe = ec3.getCodon(4, true);
        assertEquals("GGA", cwe.codon());
        assertEquals("GG", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("A", cwe.fixedSuffix());
/*
        cwe = ec.getCodon(1, false);
        assertEquals("GGG", cwe.codon());
        assertEquals("GGG", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix());

        cwe = ec2.getCodon(0, false);
        assertEquals("CCG", cwe.codon());
        assertEquals("CCG", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix());

        cwe = ec2.getCodon(4, false);
        assertEquals("TAA", cwe.codon());
        assertEquals("AA", cwe.variablePart());
        assertEquals("", cwe.fixedPrefix());
        assertEquals("T", cwe.fixedSuffix());

        cwe = ec3.getCodon(0, false);
        assertEquals("TCC", cwe.codon());
        assertEquals("CC", cwe.variablePart());
        assertEquals("T", cwe.fixedPrefix());
        assertEquals("", cwe.fixedSuffix()); */
    }

    @Test
    public void baseSequenceWithDeletionAppliedTest()
    {
        assertEquals("AAACCCGGG", ec.baseSequenceWithDeletionApplied(0, 3, true));
        assertEquals("TTTCCCGGG", ec.baseSequenceWithDeletionApplied(3, 6, true));
        assertEquals("TTACCCGGG", ec.baseSequenceWithDeletionApplied(2, 5, true));

        assertEquals("AAAACCCGG", ec2.baseSequenceWithDeletionApplied(0, 3, true));
        assertEquals("TATTAGGAT", ec6.baseSequenceWithDeletionApplied(3, 6, true));

        assertEquals("GGGTTTAAA", ec.baseSequenceWithDeletionApplied(0, 3, false));
        assertEquals("CCCTTTAAA", ec.baseSequenceWithDeletionApplied(3, 6, false));
        assertEquals("CCGTTTAAA", ec.baseSequenceWithDeletionApplied(2, 5, false));

        assertEquals("GGTTTAAAT", ec2.baseSequenceWithDeletionApplied(0, 3, false));
        assertEquals("ATCCGAATA", ec6.baseSequenceWithDeletionApplied(3, 6, false));
    }

    @Test
    public void baseSequenceWithDuplicationAppliedTest()
    {
        assertEquals("TTTTTTAAACCCGGG", ec.baseSequenceWithDuplicationApplied(0, 3, true));
        assertEquals("TTTAAAAAACCCGGG", ec.baseSequenceWithDuplicationApplied(3, 6, true));
        assertEquals("TTTAATAAACCCGGG", ec.baseSequenceWithDuplicationApplied(2, 5, true));

        assertEquals("ATTTTTTAAACCCGG", ec2.baseSequenceWithDuplicationApplied(0, 3, true));
        assertEquals("TATTAACCACCGGAT", ec6.baseSequenceWithDuplicationApplied(3, 6, true));

        assertEquals("CCCCCCGGGTTTAAA", ec.baseSequenceWithDuplicationApplied(0, 3, false));
        assertEquals("CCCGGGGGGTTTAAA", ec.baseSequenceWithDuplicationApplied(3, 6, false));
        assertEquals("CCCGGCGGGTTTAAA", ec.baseSequenceWithDuplicationApplied(2, 5, false));

        assertEquals("CCGCCGGGTTTAAAT", ec2.baseSequenceWithDuplicationApplied(0, 3, false));
        assertEquals("ATCCGGTTGTTAATA", ec6.baseSequenceWithDuplicationApplied(3, 6, false));
    }

    @Test
    public void baseSequenceWithInsertionAppliedTest()
    {
        assertEquals("GTGTTTAAACCCGGG", ec.baseSequenceWithInsertionApplied(0, "GTG", true));
        assertEquals("TGTGTTAAACCCGGG", ec.baseSequenceWithInsertionApplied(1, "GTG", true));
        assertEquals("TTTGTGAAACCCGGG", ec.baseSequenceWithInsertionApplied(3, "GTG", true));
        assertEquals("TTTAAACCCGGTACG", ec.baseSequenceWithInsertionApplied(11, "TAC", true));
        assertEquals("TTTAAACCCGGGTAC", ec.baseSequenceWithInsertionApplied(12, "TAC", true));
        assertEquals("AGGGTTTAAACCCGG", ec2.baseSequenceWithInsertionApplied(0, "GGG", true));

        assertEquals("GTACCCGGGTTTAAA", ec.baseSequenceWithInsertionApplied(0, "TAC", false));
        assertEquals("CCCGTAGGGTTTAAA", ec.baseSequenceWithInsertionApplied(3, "TAC", false));
        assertEquals("CCGAAAGGTTTAAAT", ec2.baseSequenceWithInsertionApplied(3, "TTT", false));
        assertEquals("ATCCGGTACGTAATA", ec6.baseSequenceWithInsertionApplied(5, "CGT", false));
    }

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
    public void getSplitSequenceTest()
    {
        assertEquals(new SplitCodonSequence("GGG", "", 109), ec.getSplitSequenceForCodons(4,1, true));
        assertEquals(new SplitCodonSequence("CCC", "", 106), ec.getSplitSequenceForCodons(3,1, true));
        assertEquals(new SplitCodonSequence("CCCGGG", "", 106), ec.getSplitSequenceForCodons(3,2, true));
        assertEquals(new SplitCodonSequence("TTTAAA", "", 100), ec.getSplitSequenceForCodons(1,2, true));
        assertEquals(new SplitCodonSequence("A", "TTTAA", 100), ec2.getSplitSequenceForCodons(0,2, true));
        assertEquals(new SplitCodonSequence("TAAACC", "", 102), ec2.getSplitSequenceForCodons(1,2, true));
        assertEquals(new SplitCodonSequence("TA", "TTAACCG", 100), ec6.getSplitSequenceForCodons(0,3, true));
    }

    @Test
    public void splitAtEndTest()
    {
        assertEquals(new SplitCodonSequence("GG", "A", 109), ec3.getSplitSequenceForCodons(4,1, true));
        assertEquals(new SplitCodonSequence("CCGGG", "A", 106), ec3.getSplitSequenceForCodons(3,2, true));
    }
}
