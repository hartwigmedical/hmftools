package com.hartwig.hmftools.common.sequence;

import static junit.framework.TestCase.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class SequenceAlignerTest
{
    @Before
    public void setUp()
    {
        // org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE);
    }

    @Test
    public void testAlignSequenceSimple()
    {
        String seq;
        String refSeq;

        testAlignSequenceHelper("TTAGGGTTAGGGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMMMMMMMMMMMMMMMMM");
        testAlignSequenceHelper("TTAGGTTAGCGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMM-MMMMMMSMMMMMMM");

        seq =         "TTAGGTTAGCGTTAGTG";
        refSeq = "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMM-MMMMMMSMMMM--M-M---M", 0, seq.length(), 5, refSeq.length());

        seq = "TTAGACGTC";
        refSeq = "TTAGTTAGACGTCCGTC";
        testAlignSequenceHelper(seq, refSeq, "MMMMM----MMMM", 0, seq.length(), 4, refSeq.length());

        seq = "TT";
        refSeq = "TTA";
        testAlignSequenceHelper(seq, refSeq, "MM", 0, seq.length(), 0, 2);

        seq = "TAG";
        refSeq = "TG";
        testAlignSequenceHelper(seq, refSeq, "M+M", 0, seq.length(), 0, refSeq.length());

        seq =    "TG";
        refSeq = "TAG";
        testAlignSequenceHelper(seq, refSeq, "M-M", 0, seq.length(), 0, refSeq.length());
        seq =   "TTAGG";
        refSeq = "TAGG";
        testAlignSequenceHelper(seq, refSeq, "MMMM", 1, seq.length(), 0, refSeq.length());

        seq =       "TTAGGG";
        refSeq = "CCCTAACCC";
        testAlignSequenceHelper(seq, refSeq, "MSMSSS", 3, seq.length(), 0, refSeq.length());
    }

    @Test
    public void testAlignSubsequenceSimple()
    {
        String seq;
        String refSeq;
        seq = "TTAGGGTTAGGGTTAGGG";
        refSeq = "TTAGGGTTAGGGTTAGGG";
        testAlignSubsequenceHelper(seq, refSeq, "MMMMMMMMMMMMMMMMMM", 0, seq.length(), 0, refSeq.length());
        seq = "TTAGGTTAGCGTTAGGG";
        refSeq = "TTAGGGTTAGGGTTAGGG";
        testAlignSubsequenceHelper(seq, refSeq, "MMM-MMMMMMSMMMMMMM", 0, seq.length(), 0, refSeq.length());

        seq =    "TTAGGTTAGCGTTAGTG";
        refSeq = "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG";
        testAlignSubsequenceHelper(seq, refSeq, "MMMMMMMMMSMMMMMSM", 0, seq.length(), 0, seq.length());

        // only match middle
        seq =        "TTAGACGTC";
        refSeq = "TTAGTTAGACGTCCGTC";
        testAlignSubsequenceHelper(seq, refSeq, "MMMMMMMMM", 0, seq.length(), 4, refSeq.length() - 4);

        seq =    "TT";
        refSeq = "TTA";
        testAlignSubsequenceHelper(seq, refSeq, "MM", 0, seq.length(), 0, 2);

        seq =    "TAG";
        refSeq = "TG";
        testAlignSubsequenceHelper(seq, refSeq, "M+M", 0, seq.length(), 0, refSeq.length());

        seq =    "TG";
        refSeq = "TAG";
        testAlignSubsequenceHelper(seq, refSeq, "M-M", 0, seq.length(), 0, refSeq.length());
        seq =   "TTAGG";
        refSeq = "TAGG";
        testAlignSubsequenceHelper(seq, refSeq, "MMMM", 1, seq.length(), 0, refSeq.length());

        seq =       "TTAGGG";
        refSeq = "CCCTAACCC";
        testAlignSubsequenceHelper(seq, refSeq, "MSMSSS", 0, seq.length(), 3, refSeq.length());
    }

    @Test
    public void testAlignOverlapSimple()
    {
        String leftSeq;
        String rightSeq;
        leftSeq = "ATTCCAAGAATGCGCTGTATCTGCAA";
        rightSeq =          "TGCGCTGTATCTGCAAATGGACAGTTT";
        testAlignSubsequenceHelper(leftSeq, rightSeq, "MMMMMMMMMMMMMMMM", 10, leftSeq.length(), 0, 16);
    }

    private void testAlignSequenceHelper(String seq, String refSeq, String expectedAlignOps)
    {
        testAlignSequenceHelper(seq, refSeq, expectedAlignOps, 0, seq.length(), 0, refSeq.length());
    }

    private void testAlignSequenceHelper(String seq, String refSeq, String expectedAlignOps, int seqAlignStart, int seqAlignEnd, int refSeqAlignStart, int refSeqAlignEnd)
    {
        SequenceAligner.Alignment alignment = SequenceAligner.alignSequence(seq, refSeq);
        assertEquals(expectedAlignOps, alignment.getAlignOpsString());
    }

    private void testAlignSubsequenceHelper(String seq, String refSeq, String expectedAlignOps, int seqAlignStart, int seqAlignEnd, int refSeqAlignStart, int refSeqAlignEnd)
    {
        SequenceAligner.Alignment alignment = SequenceAligner.alignSubsequence(seq, refSeq);
        assertEquals(expectedAlignOps, alignment.getAlignOpsString());
    }

    private void testAlignOverlapHelper(String leftSeq, String rightSeq, String expectedAlignOps, int seqAlignStart, int seqAlignEnd, int refSeqAlignStart, int refSeqAlignEnd)
    {
        SequenceAligner.Alignment alignment = SequenceAligner.alignOverlap(leftSeq, rightSeq);
        assertEquals(expectedAlignOps, alignment.getAlignOpsString());
    }
}