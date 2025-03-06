package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SplitCodonSequenceTest extends ReversePaveTestBase
{
    @Test
    public void couldBeDeletionInsertionTest()
    {
        assertTrue(seq("TGC", "").couldBeDeletionInsertion());
        assertTrue(seq("TGCAAT", "").couldBeDeletionInsertion());
        assertTrue(seq("TGCAA", "T").couldBeDeletionInsertion());
        assertTrue(seq("TGCA", "AT").couldBeDeletionInsertion());
        assertFalse(seq("TGC", "AAT").couldBeDeletionInsertion());
        assertTrue(seq("TG", "CAAT").couldBeDeletionInsertion());
        assertTrue(seq("T", "GCAAT").couldBeDeletionInsertion());
    }

    @Test
    public void retainedPrefixTest()
    {
        assertEquals("", seq("TTTAAACCC", "").retainedPrefix());
        assertEquals("", seq("TTTAAACC", "C").retainedPrefix());
        assertEquals("", seq("TTTAAAC", "CC").retainedPrefix());
        assertEquals("T", seq("T", "TTAAACCC").retainedPrefix());
        assertEquals("TT", seq("TT", "TAAACCC").retainedPrefix());
    }

    @Test
    public void retainedSuffixTest()
    {
        assertEquals("",  seq("TTTAAACCC", "").retainedSuffix());
        assertEquals("C", seq("TTTAAACC", "C").retainedSuffix());
        assertEquals("CC", seq("TTTAAAC", "CC").retainedSuffix());
        assertEquals("", seq("T", "TTAAACCC").retainedSuffix());
        assertEquals("", seq("TT", "TAAACCC").retainedSuffix());
    }

    @Test
    public void segmentThatIsModifiedTest()
    {
        assertEquals("TTTAAACCC", seq("TTTAAACCC", "").segmentThatIsModified());
        assertEquals("TTTAAACC", seq("TTTAAACC", "C").segmentThatIsModified());
        assertEquals("TTTAAAC", seq("TTTAAAC", "CC").segmentThatIsModified());
        assertEquals("TAAACCC", seq("TT", "TAAACCC").segmentThatIsModified());
        assertEquals("TTAAACCC", seq("T", "TTAAACCC").segmentThatIsModified());
    }

    @Test
    public void completeSequenceTest()
    {
        assertEquals("TTTAAACCC", seq("TTTAAACCC", "").completeSequence());
        assertEquals("TTTAAACCC", seq("TTTAAA", "CCC").completeSequence());
        assertEquals("TTTAAACCC", seq("TTTAA", "ACCC").completeSequence());
    }

    @Test
    public void spansTwoExonsTest()
    {
        assertFalse(seq("TTTAAACCC", "").spansTwoExons());
        assertTrue(seq("TTTAAA", "CCC").spansTwoExons());
    }
}
