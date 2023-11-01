package com.hartwig.hmftools.common.aligner;

import static junit.framework.TestCase.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class GlobalSequenceAlignerTest
{
    @Before
    public void setUp()
    {
        // org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE);
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
        testAlignSequenceHelper(seq, refSeq, "-----MMM-MMMMMMSMMMM--M-M---M");

        seq = "TTAGACGTC";
        refSeq = "TTAGTTAGACGTCCGTC";
        testAlignSequenceHelper(seq, refSeq, "----MMMMM----MMMM");

        seq = "TT";
        refSeq = "TTA";
        testAlignSequenceHelper(seq, refSeq, "MM-");

        seq = "TAG";
        refSeq = "TG";
        testAlignSequenceHelper(seq, refSeq, "M+M");

        seq =    "TG";
        refSeq = "TAG";
        testAlignSequenceHelper(seq, refSeq, "M-M");
        seq =   "TTAGG";
        refSeq = "TAGG";
        testAlignSequenceHelper(seq, refSeq, "+MMMM");

        seq =       "TTAGGG";
        refSeq = "CCCTAACCC";
        testAlignSequenceHelper(seq, refSeq, "---MSMSSS");

        seq =    "TTAGGGATT";
        refSeq = "TTAGGG";
        testAlignSequenceHelper(seq, refSeq, "MMMMMM+++");
    }

    private void testAlignSequenceHelper(String seq, String refSeq, String expectedAlignOps)
    {
        var aligner = new GlobalSequenceAligner(1, -1, -1, -1);
        GlobalSequenceAligner.Alignment alignment = aligner.alignSequence(seq, refSeq);
        assertEquals(expectedAlignOps, alignment.getOperatorsString());
    }
}