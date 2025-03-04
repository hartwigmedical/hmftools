package com.hartwig.hmftools.pave.reverse;

import static org.junit.Assert.assertEquals;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

public class PaddedExonTest extends ReversePaveTestBase
{
    PaddedExon ec1 = new PaddedExon(0,"", "", "TTTAAACCCGGG", 100, "CATG", "TACG");
    PaddedExon ec2 = new PaddedExon(2, "A", "", "TTTAAACCCGG", 100, "CATG", "TACG");
    PaddedExon ec3 = new PaddedExon(3, "", "A", "TTAAACCCGGG", 100, "CATG", "TACG");
    PaddedExon ec6 = new PaddedExon(6, "TA", "AT", "TTAACCGG", 100, "CATG", "TACG");
    PaddedExon ec7 = new PaddedExon(0,"", "", "ACGTAC", 100, "CATG", "TAGC");
    PaddedExon ec8 = new PaddedExon(0,"T", "AG", "ACGTAC", 100, "CCCC", "TTTT");

    @Test
    public void baseSequenceWithBasesReplacedAtStrandLocTest()
    {
        assertEquals("GTGCACAAACCCGGG", ec1.baseSequenceWithBasesReplacedAtStrandLocation(100,"TTT", "GTGCAC"));
        assertEquals("TGAACCCGGG", ec1.baseSequenceWithBasesReplacedAtStrandLocation(101, "TTA", "G"));
        assertEquals("TTTCCCTTTAAACCCGGG", ec1.baseSequenceWithBasesReplacedAtStrandLocation(103, "AAA", "CCCTTTAAA"));
        assertEquals("TTTAAACCCTACTAC", ec1.baseSequenceWithBasesReplacedAtStrandLocation(109, "GGG", "TACTAC"));
        assertEquals("AGGGAAACCCGG", ec2.baseSequenceWithBasesReplacedAtStrandLocation(100, "TTT", "GGG"));
        assertEquals("TATTAACTTTAAGAT", ec6.baseSequenceWithBasesReplacedAtStrandLocation(105, "CG", "TTTAA"));
    }

    @Test
    public void baseSequenceWithSingleBaseRemovedTest()
    {
        assertEquals("ACGTAT", ec7.baseSequenceWithSingleBaseRemoved(5, true));
        assertEquals("ACGTCT", ec7.baseSequenceWithSingleBaseRemoved(4, true));
        assertEquals("AGTACT", ec7.baseSequenceWithSingleBaseRemoved(1, true));
        assertEquals("CGTACT", ec7.baseSequenceWithSingleBaseRemoved(0, true));
        assertEquals("TACGACTAG", ec8.baseSequenceWithSingleBaseRemoved(3, true));

        assertEquals("GTACGC", ec7.baseSequenceWithSingleBaseRemoved(5, false));
        assertEquals("GTACTC", ec7.baseSequenceWithSingleBaseRemoved(4, false));
        assertEquals("GACGTC", ec7.baseSequenceWithSingleBaseRemoved(1, false));
        assertEquals("TACGTC", ec7.baseSequenceWithSingleBaseRemoved(0, false));
        // T  C AC TAC AG >> CTGTAGTGA
        assertEquals("CTGTAGTGA", ec8.baseSequenceWithSingleBaseRemoved(3, false));

    }

   @Test
    public void forwardStrandBaseAndLeftNeighbourTest()
    {
        assertEquals(Pair.of("G", "A"), ec7.forwardStrandBaseAndLeftNeighbour(0, true));
        assertEquals(Pair.of("A", "C"), ec7.forwardStrandBaseAndLeftNeighbour(1, true));
        assertEquals(Pair.of("A", "C"), ec7.forwardStrandBaseAndLeftNeighbour(5, true));

        // TG ACGTAC TA
        assertEquals(Pair.of("A", "C"), ec7.forwardStrandBaseAndLeftNeighbour(0, false));
        assertEquals(Pair.of("T", "A"), ec7.forwardStrandBaseAndLeftNeighbour(1, false));
        assertEquals(Pair.of("G", "A"), ec7.forwardStrandBaseAndLeftNeighbour(5, false));
    }

    @Test
    public void fromStrandCoordinatesTest()
    {
        assertEquals(0, ec1.fromStrandCoordinates(100, true));
        assertEquals(10, ec1.fromStrandCoordinates(110, true));
        assertEquals(5, ec2.fromStrandCoordinates(105, true));

        assertEquals(12, ec1.fromStrandCoordinates(100, false));
        assertEquals(2, ec1.fromStrandCoordinates(110, false));
        assertEquals(6, ec2.fromStrandCoordinates(105, false));
    }

    @Test
    public void toStrandCoordinatesTest()
    {
        assertEquals(100, ec1.toStrandCoordinates(0, true));
        assertEquals(110, ec1.toStrandCoordinates(10, true));
        assertEquals(105, ec2.toStrandCoordinates(5, true));

        assertEquals(112, ec1.toStrandCoordinates(0, false));
        assertEquals(111, ec1.toStrandCoordinates(1, false));
        assertEquals(102, ec1.toStrandCoordinates(10, false));
        assertEquals(106, ec2.toStrandCoordinates(5, false));

        for(int i=0; i<12; i++)
        {
            assertEquals(i, ec1.fromStrandCoordinates(ec1.toStrandCoordinates(i, true),true));
            assertEquals(i, ec3.fromStrandCoordinates(ec3.toStrandCoordinates(i, true),true));
            assertEquals(i, ec1.fromStrandCoordinates(ec1.toStrandCoordinates(i, false),false));
            assertEquals(i, ec3.fromStrandCoordinates(ec3.toStrandCoordinates(i, false),false));
        }
    }

    @Test
    public void baseAtTest()
    {
        assertEquals("T", ec1.baseAt(0, true));
        assertEquals("A", ec1.baseAt(3, true));
        assertEquals("G", ec1.baseAt(11, true));
        assertEquals("T", ec2.baseAt(0, true));
        assertEquals("T", ec2.baseAt(1, true));
        assertEquals("T", ec6.baseAt(0, true));
        assertEquals("A", ec6.baseAt(3, true));

        assertEquals("G", ec1.baseAt(0, false));
        assertEquals("C", ec1.baseAt(3, false));
        assertEquals("T", ec1.baseAt(11, false));
        assertEquals("G", ec2.baseAt(0, false));
        assertEquals("G", ec2.baseAt(1, false));
        assertEquals("G", ec6.baseAt(0, false));
        assertEquals("C", ec6.baseAt(3, false));
    }

    @Test
    public void codonLocationTest()
    {
        assertEquals(0, ec1.codonLocationInExonBody(1, true));
        assertEquals(3, ec1.codonLocationInExonBody(2, true));
        assertEquals(6, ec1.codonLocationInExonBody(3, true));
        assertEquals(-1, ec2.codonLocationInExonBody(0, true));
        assertEquals(2, ec2.codonLocationInExonBody(1, true));
        assertEquals(5, ec2.codonLocationInExonBody(2, true));
        assertEquals(-2, ec6.codonLocationInExonBody(0, true));
        assertEquals(1, ec6.codonLocationInExonBody(1, true));
        assertEquals(4, ec6.codonLocationInExonBody(2, true));

        // TTTAAACCCGGG
        assertEquals(9, ec1.codonLocationInExonBody(1, false));
        assertEquals(6, ec1.codonLocationInExonBody(2, false));
        assertEquals(3, ec1.codonLocationInExonBody(3, false));
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
        assertEquals("T", ec1.baseImmediatelyBefore(1));
        assertEquals("T", ec1.baseImmediatelyBefore(3));
        assertEquals("G", ec1.baseImmediatelyBefore(11));

        assertEquals("G", ec2.baseImmediatelyBefore(0));
        assertEquals("T", ec2.baseImmediatelyBefore(1));

        assertEquals("G", ec6.baseImmediatelyBefore(0));
        assertEquals("T", ec6.baseImmediatelyBefore(1));
    }

    @Test
    public void getCodonTest()
    {
        CodonWithinExons cwe = ec1.getCodon(2, true);
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
    public void baseSequenceWithFramePreservingDeletionAppliedTest()
    {
        assertEquals("AAACCCGGG", ec1.baseSequenceWithFramePreservingDeletionApplied(0, 3, true));
        assertEquals("TTTCCCGGG", ec1.baseSequenceWithFramePreservingDeletionApplied(3, 6, true));
        assertEquals("TTACCCGGG", ec1.baseSequenceWithFramePreservingDeletionApplied(2, 5, true));

        assertEquals("AAAACCCGG", ec2.baseSequenceWithFramePreservingDeletionApplied(0, 3, true));
        assertEquals("TATTAGGAT", ec6.baseSequenceWithFramePreservingDeletionApplied(3, 6, true));

        assertEquals("GGGTTTAAA", ec1.baseSequenceWithFramePreservingDeletionApplied(0, 3, false));
        assertEquals("CCCTTTAAA", ec1.baseSequenceWithFramePreservingDeletionApplied(3, 6, false));
        assertEquals("CCGTTTAAA", ec1.baseSequenceWithFramePreservingDeletionApplied(2, 5, false));

        assertEquals("GGTTTAAAT", ec2.baseSequenceWithFramePreservingDeletionApplied(0, 3, false));
        assertEquals("ATCCGAATA", ec6.baseSequenceWithFramePreservingDeletionApplied(3, 6, false));
    }

    @Test
    public void baseSequenceWithDuplicationAppliedTest()
    {
        assertEquals("TTTTTTAAACCCGGG", ec1.baseSequenceWithDuplicationApplied(0, 3, true));
        assertEquals("TTTAAAAAACCCGGG", ec1.baseSequenceWithDuplicationApplied(3, 6, true));
        assertEquals("TTTAATAAACCCGGG", ec1.baseSequenceWithDuplicationApplied(2, 5, true));

        assertEquals("ATTTTTTAAACCCGG", ec2.baseSequenceWithDuplicationApplied(0, 3, true));
        assertEquals("TATTAACCACCGGAT", ec6.baseSequenceWithDuplicationApplied(3, 6, true));

        assertEquals("CCCCCCGGGTTTAAA", ec1.baseSequenceWithDuplicationApplied(0, 3, false));
        assertEquals("CCCGGGGGGTTTAAA", ec1.baseSequenceWithDuplicationApplied(3, 6, false));
        assertEquals("CCCGGCGGGTTTAAA", ec1.baseSequenceWithDuplicationApplied(2, 5, false));

        assertEquals("CCGCCGGGTTTAAAT", ec2.baseSequenceWithDuplicationApplied(0, 3, false));
        assertEquals("ATCCGGTTGTTAATA", ec6.baseSequenceWithDuplicationApplied(3, 6, false));
    }

    @Test
    public void baseSequenceWithInsertionAppliedTest()
    {
        assertEquals("GTGTTTAAACCCGGG", ec1.baseSequenceWithInsertionApplied(0, "GTG", true));
        assertEquals("TGTGTTAAACCCGGG", ec1.baseSequenceWithInsertionApplied(1, "GTG", true));
        assertEquals("TTTGTGAAACCCGGG", ec1.baseSequenceWithInsertionApplied(3, "GTG", true));
        assertEquals("TTTAAACCCGGTACG", ec1.baseSequenceWithInsertionApplied(11, "TAC", true));
        assertEquals("TTTAAACCCGGGTAC", ec1.baseSequenceWithInsertionApplied(12, "TAC", true));
        assertEquals("AGGGTTTAAACCCGG", ec2.baseSequenceWithInsertionApplied(0, "GGG", true));

        assertEquals("GTACCCGGGTTTAAA", ec1.baseSequenceWithInsertionApplied(0, "TAC", false));
        assertEquals("CCCGTAGGGTTTAAA", ec1.baseSequenceWithInsertionApplied(3, "TAC", false));
        assertEquals("CCGAAAGGTTTAAAT", ec2.baseSequenceWithInsertionApplied(3, "TTT", false));
        assertEquals("ATCCGGTACGTAATA", ec6.baseSequenceWithInsertionApplied(5, "CGT", false));
    }

    @Test
    public void basesBetweenTest()
    {
        assertEquals("T", ec1.basesBetween(0,0));
        assertEquals("TT", ec1.basesBetween(0,1));
        assertEquals("TTTAAACCCGGG", ec1.basesBetween(0,11));
    }

    @Test
    public void basesBefore()
    {
        assertEquals(0, ec1.numberOfBasesInPreviousExon());
        assertEquals(1, ec2.numberOfBasesInPreviousExon());
        assertEquals(2, ec6.numberOfBasesInPreviousExon());
    }

    @Test
    public void getSplitSequenceTest()
    {
        assertEquals(new SplitCodonSequence("GGG", "", 109), ec1.getSplitSequenceForCodons(4,1, true));
        assertEquals(new SplitCodonSequence("CCC", "", 106), ec1.getSplitSequenceForCodons(3,1, true));
        assertEquals(new SplitCodonSequence("CCCGGG", "", 106), ec1.getSplitSequenceForCodons(3,2, true));
        assertEquals(new SplitCodonSequence("TTTAAA", "", 100), ec1.getSplitSequenceForCodons(1,2, true));
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
