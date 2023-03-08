package com.hartwig.hmftools.common.aligner;

import static junit.framework.TestCase.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class LocalSequenceAlignerTest
{
    @Before
    public void setUp()
    {
        //org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE);
    }

    @Test
    public void testAlignSequenceSimple()
    {
        String seq;
        String refSeq;

        //testAlignSequenceHelper("TTAGGGTTAGGGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMMMMMMMMMMMMMMMMM");
        //testAlignSequenceHelper("TTAGGTTAGCGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMM-MMMMMMSMMMMMMM");

        seq =    "TTAGGTTAGCGTTAGTG";
        refSeq = "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMMMMMMMMSMMMMMSM", 0, seq.length(), 0, 17);

        seq =        "TTAGACGTC";
        refSeq = "TTAGTTAGACGTCCGTC";
        testAlignSequenceHelper(seq, refSeq, "MMMMMMMMM", 0, seq.length(), 4, 13);

        seq =    "TT";
        refSeq = "TTA";
        testAlignSequenceHelper(seq, refSeq, "MM", 0, seq.length(), 0, 2);

        seq = "TTAGTT";
        refSeq = "TTGTT";
        testAlignSequenceHelper(seq, refSeq, "MM+MMM", 0, seq.length(), 0, refSeq.length());

        seq =    "TTGTT";
        refSeq = "TTAGTT";
        testAlignSequenceHelper(seq, refSeq, "MM-MMM", 0, seq.length(), 0, refSeq.length());

        seq =   "TTAGG";
        refSeq = "TAGG";
        testAlignSequenceHelper(seq, refSeq, "MMMM", 1, seq.length(), 0, refSeq.length());

        seq =      "TTAGGG";
        refSeq = "CCCTAACCC";
        testAlignSequenceHelper(seq, refSeq, "MM", 1, 3, 3, 5);
    }

    @Test
    public void testAlignSubsequenceSimple()
    {
        String seq;
        String refSeq;
        seq = "TTAGGGTTAGGGTTAGGG";
        refSeq = "TTAGGGTTAGGGTTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMMMMMMMMMMMMMMMMM", 0, seq.length(), 0, refSeq.length());
        seq =    "TTAGGTTAGCGTTAGGG";
        refSeq = "TTAGGGTTAGGGTTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMM-MMMMMMSMMMMMMM", 0, seq.length(), 0, refSeq.length());

        seq =    "TTAGGTTAGCGTTAGTG";
        refSeq = "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMMMMMMMMSMMMMMSM", 0, seq.length(), 0, seq.length());

        // only match middle
        seq =        "TTAGACGTC";
        refSeq = "TTAGTTAGACGTCCGTC";
        testAlignSequenceHelper(seq, refSeq, "MMMMMMMMM", 0, seq.length(), 4, refSeq.length() - 4);

        seq =    "TT";
        refSeq = "TTA";
        testAlignSequenceHelper(seq, refSeq, "MM", 0, seq.length(), 0, 2);

        seq =   "TTAGG";
        refSeq = "TAGG";
        testAlignSequenceHelper(seq, refSeq, "MMMM", 1, seq.length(), 0, refSeq.length());
    }

    @Test
    public void testAlignOverlapSimple()
    {
        String leftSeq;
        String rightSeq;
        leftSeq = "ATTCCAAGAATGCGCTGTATCTGCAA";
        rightSeq =          "TGCGCTGTATCTGCAAATGGACAGTTT";
        testAlignSequenceHelper(leftSeq, rightSeq, "MMMMMMMMMMMMMMMM", 10, leftSeq.length(), 0, 16);
    }

    @Test
    public void testGapPenalty()
    {
        // test the logic of gap opening and gap extension
        String left =  "AGAC-TTAGGC".replace("-", "");
        String right =  "GACTTTAGGCCT";
        LocalSequenceAligner aligner = new LocalSequenceAligner(1, -3, -2, -1);
        LocalSequenceAligner.Alignment alignment =  aligner.alignSequence(left, right);
        assertEquals("MMM-MMMMMM", alignment.getOperatorsString());

        // check that score is correct, there are 9 matching bases and 1 gap = 9 - 2 = 7
        assertEquals(7, alignment.getScore());

        // check reverse works the same
        alignment =  aligner.alignSequence(right, left);
        assertEquals("MMM+MMMMMM", alignment.getOperatorsString());

        // check that score is correct, there are 9 matching bases and 1 gap = 9 - 2 = 7
        assertEquals(7, alignment.getScore());


        // now test gap extension
        left =  "AGACTT---AGGCC".replace("-", "");
        right =  "GACTTCTGAGGCCT";
        alignment =  aligner.alignSequence(left, right);
        assertEquals("MMMMM---MMMMM", alignment.getOperatorsString());

        // check that score is correct, there are 10 matching bases and 1 gap open and 2 gap extend = 10 - 2 - 2 = 6
        assertEquals(6, alignment.getScore());

        // check that reverse work the same
        alignment =  aligner.alignSequence(right, left);
        assertEquals("MMMMM+++MMMMM", alignment.getOperatorsString());

        // check that score is correct, there are 10 matching bases and 1 gap open and 2 gap extend = 10 - 2 - 2 = 6
        assertEquals(6, alignment.getScore());
    }

    private void testAlignSequenceHelper(String seq, String refSeq, String expectedAlignOps, int seq1AlignStart, int seq1AlignEnd,
            int seq2AlignStart, int seq2AlignEnd)
    {
        LocalSequenceAligner aligner = new LocalSequenceAligner(2, -1, -2, -1);
        LocalSequenceAligner.Alignment alignment =  aligner.alignSequence(seq, refSeq);
        assertEquals(expectedAlignOps, alignment.getOperatorsString());
        assertEquals(seq1AlignStart, alignment.getFirstSequenceAlignStart());
        assertEquals(seq1AlignEnd, alignment.getFirstSequenceAlignEnd());
        assertEquals(seq2AlignStart, alignment.getSecondSequenceAlignStart());
        assertEquals(seq2AlignEnd, alignment.getSecondSequenceAlignEnd());
    }
}