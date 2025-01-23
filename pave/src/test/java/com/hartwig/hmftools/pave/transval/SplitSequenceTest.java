package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SplitSequenceTest extends TransvalTestBase
{
    @Test
    public void couldBeDeletionInsertionTest()
    {
        assertTrue(seq("TGC", null).couldBeDeletionInsertion());
        assertTrue(seq("TGCAAT", null).couldBeDeletionInsertion());
        assertTrue(seq("TGCAA", "T").couldBeDeletionInsertion());
        assertTrue(seq("TGCA", "AT").couldBeDeletionInsertion());
        assertFalse(seq("TGC", "AAT").couldBeDeletionInsertion());
        assertTrue(seq("TG", "CAAT").couldBeDeletionInsertion());
        assertTrue(seq("T", "GCAAT").couldBeDeletionInsertion());
    }

    @Test
    public void segmentThatIsModifiedTest()
    {
        assertEquals("TTTAAACCC", seq("TTTAAACCC", null).segmentThatIsModified());
        assertEquals("TTTAAACC", seq("TTTAAACC", "C").segmentThatIsModified());
        assertEquals("TTTAAAC", seq("TTTAAAC", "CC").segmentThatIsModified());
        assertEquals("TAAACCC", seq("TT", "TAAACCC").segmentThatIsModified());
        assertEquals("TTAAACCC", seq("T", "TTAAACCC").segmentThatIsModified());
    }

    @Test
    public void completeSequenceTest()
    {
        assertEquals("TTTAAACCC", seq("TTTAAACCC", null).completeSequence());
        assertEquals("TTTAAACCC", seq("TTTAAA", "CCC").completeSequence());
        assertEquals("TTTAAACCC", seq("TTTAA", "ACCC").completeSequence());
    }

    @Test
    public void spansTwoExonsTest()
    {
        assertFalse(seq("TTTAAACCC", null).spansTwoExons());
        assertTrue(seq("TTTAAA", "CCC").spansTwoExons());
    }
}
