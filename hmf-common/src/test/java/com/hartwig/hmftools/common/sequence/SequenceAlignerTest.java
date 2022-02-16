package com.hartwig.hmftools.common.sequence;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import org.apache.logging.log4j.Level;
import org.junit.Before;
import org.junit.Test;

public class SequenceAlignerTest
{
    @Before
    public void setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE);
    }

    @Test
    public void testAlignSequenceSimple()
    {
        testAlignSequenceHelper("TTAGGGTTAGGGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMMMMMMMMMMMMMMMMM");
        testAlignSequenceHelper("TTAGGTTAGCGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMM-MMMMMMSMMMMMMM");

        testAlignSequenceHelper("TTAGGTTAGCGTTAGTG", "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG", "-----MMM-MMMMMMSMMMM--M-M---M");


        testAlignSequenceHelper("TTAGACGTC", "TTAGTTAGACGTCCGTC", "----MMMMM----MMMM");

        testAlignSequenceHelper("TT", "TTA", "MM-");

        testAlignSequenceHelper("TAG", "TG", "M+M");

        testAlignSequenceHelper("TG", "TAG", "M-M");
        testAlignSequenceHelper("TTAGG", "TAGG", "+MMMM");

        testAlignSequenceHelper("TTAGGG", "CCCTAACCC", "---MSMSSS");
    }

    @Test
    public void testAlignSubsequenceSimple()
    {
        testAlignSubsequenceHelper("TTAGGGTTAGGGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMMMMMMMMMMMMMMMMM");
        testAlignSubsequenceHelper("TTAGGTTAGCGTTAGGG", "TTAGGGTTAGGGTTAGGG", "MMM-MMMMMMSMMMMMMM");

        testAlignSubsequenceHelper("TTAGGTTAGCGTTAGTG", "TTAGGTTAGGGTTAGGGTTAGGGTTAGGG", "MMMMMMMMMSMMMMMSM------------");


        testAlignSubsequenceHelper("TTAGACGTC", "TTAGTTAGACGTCCGTC", "----MMMMMMMMM----");

        testAlignSubsequenceHelper("TT", "TTA", "MM-");

        testAlignSubsequenceHelper("TAG", "TG", "M+M");

        testAlignSubsequenceHelper("TG", "TAG", "M-M");
        testAlignSubsequenceHelper("TTAGG", "TAGG", "+MMMM");

        testAlignSubsequenceHelper("TTAGGG", "CCCTAACCC", "---MSMSSS");
    }

    private void testAlignSequenceHelper(String seq, String refSeq, String expectedAlignOps)
    {
        List<SequenceAligner.AlignOp> alignOps = SequenceAligner.alignSequence(seq, refSeq);
        assertEquals(expectedAlignOps, SequenceAligner.AlignOp.toString(alignOps));
    }

    private void testAlignSubsequenceHelper(String seq, String refSeq, String expectedAlignOps)
    {
        List<SequenceAligner.AlignOp> alignOps = SequenceAligner.alignSubsequence(seq, refSeq);
        assertEquals(expectedAlignOps, SequenceAligner.AlignOp.toString(alignOps));
    }
}

